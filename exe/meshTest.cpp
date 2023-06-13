#include <iostream>
#include "../src/ConvexHull.h"
#include "../src/Stencil.h"
using namespace std;



void ConvexHullCube()
{
    std::vector<Vector3> vertices;
    vertices.push_back(Vector3(0,0,0));
    vertices.push_back(Vector3(0,0,1));
    vertices.push_back(Vector3(0,1,0));
    vertices.push_back(Vector3(0,1,1));
    vertices.push_back(Vector3(1,0,0));
    vertices.push_back(Vector3(1,0,1));
    vertices.push_back(Vector3(1,1,0));
    vertices.push_back(Vector3(1,1,1));

    ConvexHull convexHull(vertices);
    Mesh mesh = convexHull.GetMesh();
    mesh.WriteToCsv(OUTPUTDIR,"cube");
}
void ConvexHullStencil()
{
    for(int i=3; i<=31; i+=2)
    {
        LebedevStencil stencil(i);
        std::vector<Vector3> vertices;
        vertices.reserve(stencil.nDir);

        for(int d=0; d<stencil.nDir; d++)
            vertices.push_back(stencil.Cv3(d));
        
        ConvexHull convexHull(vertices);
        Mesh mesh = convexHull.GetMesh();
        // mesh.WriteToUnityCsv("output","Lebedev" + to_string(i));

        convexHull.OriginalOrdering(vertices);
        mesh = convexHull.GetMesh();
        mesh.WriteToCsv(OUTPUTDIR,"Lebedev" + to_string(i) + "ordered");
    }
}



int main()
{
    // ConvexHullCube();
    ConvexHullStencil();
}