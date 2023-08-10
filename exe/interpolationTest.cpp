#include <iostream>
#include "../src/Grid.h"
#include "../src/Interpolation.h"
#include "../src/Stencil.h"
#include "../src/InterpolationGrid.h"
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
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        data0[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 1,0,0,0,0,0,0,1);
        data1[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,1,0,0,0,0,1,0);
        data2[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,0,1,0,0,1,0,0);
        data3[ijk] = 1;
    }
    grid.WriteFrametoCsv(0,data0,data1,data2,data3,OUTPUTDIR,"TrilinearInterpolation");
    
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



void TestStencilInterpolations(Stencil stencilSrc, InterpolationGrid interpGrid)
{
    GaussLegendreStencil stencilDst(41);
    double* valueSrc = new double[stencilSrc.nDir]();
    double* valueDst1 = new double[stencilDst.nDir]();
    double* valueDst2 = new double[stencilDst.nDir]();
    double* valueDst3 = new double[stencilDst.nDir]();
    double* valueDst4 = new double[stencilDst.nDir]();

    // Fill source array:
    for(size_t d=0; d<stencilSrc.nDir; d++)
        valueSrc[d] = RandomRange(1.0,1.3);

	Profiler::Session& session = Profiler::Session::Get();
    // Barycentric Interpolate:
    {
        PROFILE_SCOPE("Barycentric Interpolation");
        for(size_t d=0; d<stencilDst.nDir; d++)
        {
            Vector3 p = stencilDst.Cv3(d);
            std::tuple<Vector3Int,Vector3> neighboursAndWeights = stencilSrc.BarycentricNeighboursAndWeights(p);
            Vector3Int neighbours = std::get<0>(neighboursAndWeights);
            Vector3 weights = std::get<1>(neighboursAndWeights);
            valueDst1[d] += weights[0] * valueSrc[neighbours[0]] + weights[1] * valueSrc[neighbours[1]] + weights[2] * valueSrc[neighbours[2]];
        }
    }

    // Voronoi Interpolate:
    {
        PROFILE_SCOPE("Voronoi Interpolate");
        for(size_t d=0; d<stencilDst.nDir; d++)
        {
            Vector3 p = stencilDst.Cv3(d);
            std::tuple<std::vector<size_t>, std::vector<double>> neighboursAndWeights = stencilSrc.VoronoiNeighboursAndWeights(p);
            std::vector<size_t> neighbours = std::get<0>(neighboursAndWeights);
            std::vector<double> weights = std::get<1>(neighboursAndWeights);
            for(int i=0; i<neighbours.size(); i++)
                valueDst2[d] += weights[i] * valueSrc[neighbours[i]];
        }
    }

    // Optimized Barycentric Interpolate:
    {
        PROFILE_SCOPE("Optimized Barycentric Interpolate");
        for(size_t d=0; d<stencilDst.nDir; d++)
        {
            Tensor3 c = stencilDst.Ct3(d);
            double i = interpGrid.i(c);
            double j = interpGrid.j(c);
            size_t i0 = std::floor(i); size_t i1 = i0 + 1;
            size_t j0 = std::floor(j); size_t j1 = (j0 + 1) % interpGrid.nPh;

            Vector3Int neighbours00 = interpGrid.BarycentricNeighbours(i0,j0);
            Vector3Int neighbours01 = interpGrid.BarycentricNeighbours(i0,j1);
            Vector3Int neighbours10 = interpGrid.BarycentricNeighbours(i1,j0);
            Vector3Int neighbours11 = interpGrid.BarycentricNeighbours(i1,j1);
            Vector3 weights00 = interpGrid.BarycentricWeights(i0,j0);
            Vector3 weights01 = interpGrid.BarycentricWeights(i0,j1);
            Vector3 weights10 = interpGrid.BarycentricWeights(i1,j0);
            Vector3 weights11 = interpGrid.BarycentricWeights(i1,j1);
            double value00 = weights00[0] * valueSrc[neighbours00[0]] + weights00[1] * valueSrc[neighbours00[1]] + weights00[2] * valueSrc[neighbours00[2]];
            double value01 = weights01[0] * valueSrc[neighbours01[0]] + weights01[1] * valueSrc[neighbours01[1]] + weights01[2] * valueSrc[neighbours01[2]];
            double value10 = weights10[0] * valueSrc[neighbours10[0]] + weights10[1] * valueSrc[neighbours10[1]] + weights10[2] * valueSrc[neighbours10[2]];
            double value11 = weights11[0] * valueSrc[neighbours11[0]] + weights11[1] * valueSrc[neighbours11[1]] + weights11[2] * valueSrc[neighbours11[2]];
            valueDst3[d] = BilinearInterpolation(i-i0,j-j0,value00,value01,value10,value11);
        }
    }

    // Optimized Voronoi Interpolate:
    {
        PROFILE_SCOPE("Optimized Voronoi Interpolate");
        for(size_t d=0; d<stencilDst.nDir; d++)
        {
            Tensor3 c = stencilDst.Ct3(d);
            double i = interpGrid.i(c);
            double j = interpGrid.j(c);
            size_t i0 = std::floor(i); size_t i1 = i0 + 1;
            size_t j0 = std::floor(j); size_t j1 = (j0 + 1) % interpGrid.nPh;

            std::span<const size_t> neighbours00 = interpGrid.VoronoiNeighbours(i0,j0);
            std::span<const size_t> neighbours01 = interpGrid.VoronoiNeighbours(i0,j1);
            std::span<const size_t> neighbours10 = interpGrid.VoronoiNeighbours(i1,j0);
            std::span<const size_t> neighbours11 = interpGrid.VoronoiNeighbours(i1,j1);
            std::span<const double> weights00 = interpGrid.VoronoiWeights(i0,j0);
            std::span<const double> weights01 = interpGrid.VoronoiWeights(i0,j1);
            std::span<const double> weights10 = interpGrid.VoronoiWeights(i1,j0);
            std::span<const double> weights11 = interpGrid.VoronoiWeights(i1,j1);
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
            valueDst4[d] = BilinearInterpolation(i-i0,j-j0,value00,value01,value10,value11);
        }
    }
	session.End();
	std::vector<std::string> names = session.GetAllFunctionNames();
	for(int i=0; i<names.size(); i++)
        session.PrintFunctionDuration(names[i]);

    // Write to file:
    ofstream file0(OUTPUTDIR + "Src.txt");
    ofstream file1(OUTPUTDIR + "DstBarycentric.txt");
    ofstream file2(OUTPUTDIR + "DstVoronoi.txt");
    ofstream file3(OUTPUTDIR + "DstBarycentricFast.txt");
    ofstream file4(OUTPUTDIR + "DstVoronoiFast.txt");
    file0 << "x,y,z,color" << endl;
    file1 << "x,y,z,color" << endl;
    file2 << "x,y,z,color" << endl;
    file3 << "x,y,z,color" << endl;
    file4 << "x,y,z,color" << endl;
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
        Tensor3 c4 = valueDst4[d] * stencilDst.Ct3(d);
        file1 << c1[1] << "," << c1[2] << "," << c1[3] << "," << 0.25 << "\n";
        file2 << c2[1] << "," << c2[2] << "," << c2[3] << "," << 0.50 << "\n";
        file3 << c3[1] << "," << c3[2] << "," << c3[3] << "," << 0.75 << "\n";
        file4 << c4[1] << "," << c4[2] << "," << c4[3] << "," << 1.00 << "\n";
    }
    file0.close();
    file1.close();
    file2.close();
    file3.close();
    file4.close();
    
    cout << "Three files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



int main()
{
    // TestTrilinearInterpolation();
    // TestRayTriangleIntersection();
    // TestStencilInterpolations(LebedevStencil(21));
    LebedevStencil stencil(21,3,8,M_PI/8.0);
    InterpolationGrid interpGrid(200,500,stencil);
    TestStencilInterpolations(stencil,interpGrid);
}