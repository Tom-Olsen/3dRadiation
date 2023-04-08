#include <iostream>
#include "../src/Grid.h"
#include "../src/Interpolation.hh"
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
    
    cout << (string)OUTPUTDIR + "TrilinearInterpolation... has been created. Plot it with ParaView (Filter:Table to Points)." << endl;
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
    Tensor3 weights;
    int n = 101;

    ofstream file((string)OUTPUTDIR + (string)"BarycentricInterpolation.csv");
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
            Tensor3 color = color0 * weights[1] + color1 * weights[2] + color2 * weights[3];
            file << intersection[1] << "," << intersection[2] << "," << intersection[3] << ","
                 << color[1] << "," << color[2] << "," << color[3] << "\n";
        }
    }
    file.close();

    cout << (string)OUTPUTDIR + "BarycentricInterpolation.csv has been created. Plot it with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



void Test_SphericalBarycentricWeights()
{
    Tensor3 a(1,0,0);
    Tensor3 b(0,1,0);
    Tensor3 c(0,0,1);
    
    Tensor3 p(1,1,1);
    double norm = p.EuklNorm();
    p[1] /= norm;
    p[2] /= norm;
    p[3] /= norm;

    Tensor3 weights;
    SphericalBarycentricWeights(p,a,b,c,weights);
    cout << "The results are wrong. Spherical Barycentric weights and not used yet." << endl;
    cout << endl;
}



void BarycentricWeightsTest(int nOrder)
{
    LebedevStencil stencil(nOrder);

    // Construct mesh:
    std::vector<Vector3> vertices;
    for(int d=0; d<stencil.nDir; d++)
    {
        Tensor3 c = stencil.C(d);
        Vector3 v(c[1],c[2],c[3]);
        vertices.push_back(v);
    }
    ConvexHull convexHull(vertices);
    convexHull.OriginalOrdering(vertices);
    Mesh mesh = convexHull.GetMesh();
    std::vector<Vector3Int> triangles = mesh.GetTriangles();

    ofstream file((string)OUTPUTDIR + "BarycentricWeightsTest.txt");
    file << "x,y,z,value" << endl;
    int nPhi = 30;
    int nTheta = 60;
    for (int j = 0; j < nPhi; j++)
    for (int i = 0; i < nTheta; i++)
    {
        Vector3 origin(0,0,0);
        double theta = M_PI * (i + 0.5) / (double)nTheta;
        double phi = 2.0 * M_PI * j / (double)nPhi;
        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);
        Vector3 direction(x, y, z);
        Vector3 weights, intersectionPoint;
        for (int d = 0; d < triangles.size(); d++)
        {
            Vector3Int triangle = triangles[d];
            Vector3 v0 = vertices[triangle[0]];
            Vector3 v1 = vertices[triangle[1]];
            Vector3 v2 = vertices[triangle[2]];
            if (RayTriangleIntersection(origin,direction,v0,v1,v2,intersectionPoint))
            {
                double value = 0;
                file << intersectionPoint[0] << "," << intersectionPoint[1] << "," << intersectionPoint[2] << "," << 0 << endl;
                break;
            }
            //if (BarycentricWeights(origin, direction, v0, v1, v2, weights))
            //{
            //    RayTriangleIntersection(origin,direction,v0,v1,v2,intersectionPoint);
            //    double value = 0;
            //    file << intersectionPoint[0] << "," << intersectionPoint[1] << "," << intersectionPoint[2] << "," << 0 << endl;
            //    break;
            //}
        }
    }
    file.close();
    
    cout << "A file has been created in '" + (string)OUTPUTDIR + "'. Plot it with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



int main()
{
    TestTrilinearInterpolation();
    TestRayTriangleIntersection();
    Test_SphericalBarycentricWeights();
    BarycentricWeightsTest(3);
}