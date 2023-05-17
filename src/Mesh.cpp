#include "Mesh.h"



// Constructors:
Mesh::Mesh()
{
    vertexCount = 0;
}
Mesh::Mesh(Vector3* vertices, int vertexCount, Vector3Int* triangles, int triangleCount)
{
    this->vertexCount = vertexCount;
    this->m_vertices.reserve(vertexCount);
    for(int i=0; i<vertexCount; i++)
        this->m_vertices.push_back(vertices[i]);
        
    this->m_triangles.reserve(triangleCount);
    for(int i=0; i<triangleCount; i++)
        this->m_triangles.push_back(triangles[i]);
}
Mesh::Mesh(std::vector<Vector3> vertices, std::vector<Vector3Int> triangles)
{
    this->vertexCount = vertices.size();
    this->m_vertices = vertices;
    this->m_triangles = triangles;
}

// Getters:
std::vector<Vector3> Mesh::GetVertices()
{ return m_vertices; }
std::vector<Vector3Int> Mesh::GetTriangles()
{ return m_triangles; }

// Output:
void Mesh::WriteToCsv(std::string directory, std::string name)
{
    CreateDirectory(directory);
    std::string fullPath = (directory == "") ? name + ".csv" : directory + "/" + name + ".csv";
    
    std::ofstream fileOut(fullPath);
    fileOut << "vertices:\n";
    for(const Vector3& v : m_vertices)
        fileOut << v << "\n";
    
    fileOut << "\n";
    fileOut << "triangles:\n";
    for(const Vector3Int& t : m_triangles)
        fileOut << t << "\n";
    fileOut.close();
}