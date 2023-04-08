#ifndef __INCLUDE_GUARD_Mesh_h__
#define __INCLUDE_GUARD_Mesh_h__
#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include "Utility.hh"
#include "Vector3.h"
#include "Vector3Int.h"



struct Mesh
{
private:
    std::vector<Vector3> m_vertices;
    std::vector<Vector3Int> m_triangles;
public:
    int vertexCount;

    // Constructors:
    Mesh(Vector3* m_vertices, int vertexCount, Vector3Int* m_triangles, int triangleCount);
    Mesh(std::vector<Vector3> vertices, std::vector<Vector3Int> triangles);

    // Getters:
    std::vector<Vector3> GetVertices();
    std::vector<Vector3Int> GetTriangles();

    // Output:
    void WriteToCsv(std::string directory, std::string name);
};
#endif //__INCLUDE_GUARD_Mesh_h__