#ifndef __INCLUDE_GUARD_ConvexHull_h__
#define __INCLUDE_GUARD_ConvexHull_h__
#include <iostream>
#include <vector>
#include <set>
#include "Utility.hh"
#include "Vector3.h"
#include "Vector3Int.h"
#include "Vector2Int.h"
#include "Mesh.h"
#include "DataTypes.hh"



struct ConvexHull
{
private:
    std::vector<Vector3> m_vertices;
    std::vector<Vector3> m_hull;
    std::vector<Vector3Int> m_triangles;
    Vector3 center;

public:
    ConvexHull(const std::vector<Vector3>& vertices);
    void OriginalOrdering(const std::vector<Vector3>& vertices);

    std::vector<Vector3> GetVertices();
    std::vector<Vector3Int> GetTriangles();
    Mesh GetMesh();

    void InitialTetrahedron();
    void AddVertex(const Vector3& newVertex);
    bool CanSeeTriangle(Vector3 vertex, Vector3Int triangle);
    Vector3Int OrientedTriangle(int a, int b, int c);
    bool InsideHull(Vector3 vertex);
    
    std::vector<Vector3Int> RemoveBadTriangles(Vector3 vertex);
    std::vector<Vector2Int> PolygonEdges(const std::vector<Vector3Int>& badTriangles);
    void RemoveBadVertices(const std::vector<Vector3Int>& badTriangles, std::vector<Vector2Int>& polygonEdges);
};
#endif //__INCLUDE_GUARD_ConvexHull_h__