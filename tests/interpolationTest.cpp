#include <iostream>
#include "../src/Grid.h"
#include "../src/Interpolation.h"
#include "../src/Stencil.h"
using namespace std;



void TestTrilinearInterpolation()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    RealBuffer data0(grid.nxyz);
    RealBuffer data1(grid.nxyz);
    RealBuffer data2(grid.nxyz);
    RealBuffer data3(grid.nxyz);
    for(size_t i=0; i<grid.nx; i++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t k=0; k<grid.nz; k++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        data0[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 1,0,0,0,0,0,0,1);
        data1[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,1,0,0,0,0,1,0);
        data2[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,0,1,0,0,1,0,0);
        data3[ijk] = 1;
    }
    grid.WriteFrametoCsv(0,data0,data1,data2,data3,0,OUTPUTDIR,"TrilinearInterpolation");
    
    cout << OUTPUTDIR + "TrilinearInterpolation... has been created. Plot it with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



void TestRayTriangleIntersection()
{
    Tensor3 v0( 0, 1, 1);
    Tensor3 v1(-1,-1, 1);
    Tensor3 v2( 1,-1, 1);

    Tensor3 color0(1,0,0);
    Tensor3 color1(0,1,0);
    Tensor3 color2(0,0,1);

    Tensor3 origin(0,0,0);
    Tensor3 intersection;
    Vector3 weights;
    int n = 101;

    ofstream file(OUTPUTDIR + "BarycentricInterpolation.txt");
    file << "x,y,z,r,g,b\n";
    for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    {
        double x = 2.0 * i / (n - 1.0) - 1.0;
        double y = 2.0 * j / (n - 1.0) - 1.0;
        Tensor3 direction(x,y,1);   // does not need to be normalized.

        RayTriangleIntersection(origin, direction, v0, v1, v2, intersection);
        if(BarycentricWeights(origin, direction, v0, v1, v2, weights))
        {
            Tensor3 color = color0 * weights[0] + color1 * weights[1] + color2 * weights[2];
            file << intersection[1] << "," << intersection[2] << "," << intersection[3] << ","
                 << color[1] << "," << color[2] << "," << color[3] << "\n";
        }
    }
    file.close();

    cout << OUTPUTDIR + "BarycentricInterpolation.txt has been created. Plot it with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



void TestStencilInterpolations(int nOrder)
{
    LebedevStencil stencilSrc(nOrder);
    GaussLegendreStencil stencilDst(31);
    double* valueSrc = new double[stencilSrc.nDir]();
    double* valueDst1 = new double[stencilDst.nDir]();
    double* valueDst2 = new double[stencilDst.nDir]();
    double* valueDst3 = new double[stencilDst.nDir]();

    // Fill source array:
    for(size_t d=0; d<stencilSrc.nDir; d++)
        valueSrc[d] = RandomRange(1.0,1.3);

    // Barycentric interpolate:
    for(size_t d=0; d<stencilDst.nDir; d++)
    {
        Tensor3 c = stencilDst.Ct3(d);
        Vector3Int triangle;
        Vector3 weights = stencilSrc.BarycentricWeights(c,triangle);
        valueDst1[d] = weights[0] * valueSrc[triangle[0]] + weights[1] * valueSrc[triangle[1]] + weights[2] * valueSrc[triangle[2]];
    }

    // Voronoi interpolate:
    for(size_t d=0; d<stencilDst.nDir; d++)
    {
        Tensor3 c = stencilDst.Ct3(d);
        std::vector<size_t> neighbours;
        std::vector<double> weights = stencilSrc.VoronoiWeights(c,neighbours);
        for(size_t i=0; i<weights.size(); i++)
            valueDst2[d] += weights[i] * valueSrc[neighbours[i]];
    }

    // Optimized Voronoi interpolate:
    for(size_t d=0; d<stencilDst.nDir; d++)
    {
        Tensor3 c = stencilDst.Ct3(d);
        double i = stencilSrc.sphereGrid.i(c.Theta());
        double j = std::fmod((stencilSrc.sphereGrid.j(c.Phi()) + stencilSrc.sphereGrid.nPh), stencilSrc.sphereGrid.nPh);
        size_t i0 = std::floor(i); size_t i1 = i0 + 1;
        size_t j0 = std::floor(j); size_t j1 = (j0 + 1) % stencilSrc.sphereGrid.nPh;

        std::span<const size_t> neighbours00 = stencilSrc.voronoiNeighboursOnGrid.Row(stencilSrc.sphereGrid.Index(i0,j0));
        std::span<const size_t> neighbours01 = stencilSrc.voronoiNeighboursOnGrid.Row(stencilSrc.sphereGrid.Index(i0,j1));
        std::span<const size_t> neighbours10 = stencilSrc.voronoiNeighboursOnGrid.Row(stencilSrc.sphereGrid.Index(i1,j0));
        std::span<const size_t> neighbours11 = stencilSrc.voronoiNeighboursOnGrid.Row(stencilSrc.sphereGrid.Index(i1,j1));
        std::span<const double> weights00 = stencilSrc.voronoiWeightsOnGrid.Row(stencilSrc.sphereGrid.Index(i0,j0));
        std::span<const double> weights01 = stencilSrc.voronoiWeightsOnGrid.Row(stencilSrc.sphereGrid.Index(i0,j1));
        std::span<const double> weights10 = stencilSrc.voronoiWeightsOnGrid.Row(stencilSrc.sphereGrid.Index(i1,j0));
        std::span<const double> weights11 = stencilSrc.voronoiWeightsOnGrid.Row(stencilSrc.sphereGrid.Index(i1,j1));
        double value00 = 0;
        double value01 = 0;
        double value10 = 0;
        double value11 = 0;
        for(size_t i=0; i<weights00.size(); i++)
            value00 += weights00[i] * valueSrc[neighbours00[i]];
        for(size_t i=0; i<weights01.size(); i++)
            value01 += weights01[i] * valueSrc[neighbours01[i]];
        for(size_t i=0; i<weights10.size(); i++)
            value10 += weights10[i] * valueSrc[neighbours10[i]];
        for(size_t i=0; i<weights11.size(); i++)
            value11 += weights11[i] * valueSrc[neighbours11[i]];
        valueDst3[d] = BilinearInterpolation(i-i0,j-j0,value00,value01,value10,value11);
    }

    // Write to file:
    ofstream file0(OUTPUTDIR + "Src.txt");
    ofstream file1(OUTPUTDIR + "DstBarycentric.txt");
    ofstream file2(OUTPUTDIR + "DstVoronoi.txt");
    ofstream file3(OUTPUTDIR + "DstVoronoi2.txt");
    file0 << "x,y,z,color" << endl;
    file1 << "x,y,z,color" << endl;
    file2 << "x,y,z,color" << endl;
    file3 << "x,y,z,color" << endl;
    for(int d=0; d<stencilSrc.nDir; d++)
    {
        Tensor3 c0 = valueSrc[d] * stencilSrc.Ct3(d);
        file0 << c0[1] << "," << c0[2] << "," << c0[3] << "," << 0.00 << "\n";
    }
    for(int d=0; d<stencilDst.nDir; d++)
    {
        Tensor3 c1 = valueDst1[d] * stencilDst.Ct3(d);
        Tensor3 c2 = valueDst2[d] * stencilDst.Ct3(d);
        Tensor3 c3 = valueDst3[d] * stencilDst.Ct3(d);
        file1 << c1[1] << "," << c1[2] << "," << c1[3] << "," << 0.33 << "\n";
        file2 << c2[1] << "," << c2[2] << "," << c2[3] << "," << 0.66 << "\n";
        file3 << c3[1] << "," << c3[2] << "," << c3[3] << "," << 1.00 << "\n";
    }
    file0.close();
    file1.close();
    file2.close();
    file3.close();
    
    cout << "Three files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



int main()
{
    // TestTrilinearInterpolation();
    // TestRayTriangleIntersection();
    TestStencilInterpolations(35);
}