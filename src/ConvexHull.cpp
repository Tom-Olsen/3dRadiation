#include "ConvexHull.h"



ConvexHull::ConvexHull(const std::vector<Vector3>& vertices) : m_vertices(vertices)
{
    if(m_vertices.size() < 4)
        ExitOnError("At least 4 points reqired.");
    m_hull.reserve(m_vertices.size());
    InitialTetrahedron();
    
    for(auto vertex=m_vertices.rbegin(); vertex<m_vertices.rend(); vertex++)
    {
        if (!InsideHull(*vertex))
            AddVertex(*vertex);
        m_vertices.pop_back();
    }
    m_hull.shrink_to_fit();
}

void ConvexHull::OriginalOrdering(const std::vector<Vector3>& vertices)
{
    if(vertices.size() < m_hull.size())
        ExitOnError("Reordering makes no sense, given vertices are to short.");

    // Only concider verticies that are on the convex hull:
    std::vector<Vector3> orderedHull;
    orderedHull.reserve(m_hull.size());
    std::vector<Vector2Int> permutation;
    permutation.reserve(m_hull.size());

    int newIndex = 0;
    for(int i=0; i<vertices.size(); i++)
    {
        Vector3 v = vertices[i];
        for(int j=0; j<m_hull.size(); j++)
        {
            Vector3 w = m_hull[j];
            if(v == w)
            {
                orderedHull.push_back(v);
                permutation.push_back(Vector2Int(j,newIndex));
                newIndex++;
                break;
            }
        }
    }

    // Order triangles such that the same triangles are drawn with the new ordering:
    std::vector<Vector3Int> orderedTriangles;
    orderedTriangles.reserve(m_triangles.size());
    for(const Vector3Int& triangle : m_triangles)
    {
        int a = triangle[0];
        int b = triangle[1];
        int c = triangle[2];
        for(Vector2Int indexes : permutation)
            if(a==indexes[0])
            {
                a = indexes[1];
                break;
            }
        for(Vector2Int indexes : permutation)
            if(b==indexes[0])
            {
                b = indexes[1];
                break;
            }
        for(Vector2Int indexes : permutation)
            if(c==indexes[0])
            {
                c = indexes[1];
                break;
            }
        orderedTriangles.push_back(Vector3Int(a,b,c));
    }
    m_hull = orderedHull;
    m_triangles = orderedTriangles;
}

std::vector<Vector3> ConvexHull::GetVertices()
{ return m_hull; }
std::vector<Vector3Int> ConvexHull::GetTriangles()
{ return m_triangles; }
Mesh ConvexHull::GetMesh()
{ return Mesh(m_hull, m_triangles); }

void ConvexHull::InitialTetrahedron()
{
    // Search fourth vertex, such that vertex 0,1,2,i are not coplanar:
    int vertex0 = m_vertices.size()-1;
    int vertex1 = m_vertices.size()-2;
    int vertex2 = m_vertices.size()-3;
    int vertex3 = -1;
    for(int i=m_vertices.size()-4; i>-1; i--)
    {
        if (!Vector3::AreCoplanar(m_vertices[vertex0], m_vertices[vertex1], m_vertices[vertex2], m_vertices[i]))
        {
            vertex3 = i;
            break;
        }
    }
    if (vertex3 == -1)
        ExitOnError("All given points are coplanar.");
    // Add the 4 none coplanar points to Convex m_hull to form a Tetrahedron:
    m_hull.push_back(m_vertices[vertex0]);
    m_hull.push_back(m_vertices[vertex1]);
    m_hull.push_back(m_vertices[vertex2]);
    m_hull.push_back(m_vertices[vertex3]);
    center = Vector3::GetCenter(m_hull);
    // Add Faces to Tetrahedron:
    m_triangles.push_back(OrientedTriangle(0, 1, 2));
    m_triangles.push_back(OrientedTriangle(0, 1, 3));
    m_triangles.push_back(OrientedTriangle(0, 2, 3));
    m_triangles.push_back(OrientedTriangle(1, 2, 3));
    // Remove already used m_vertices:
    m_vertices.erase(m_vertices.begin() + vertex0);
    m_vertices.erase(m_vertices.begin() + vertex1);
    m_vertices.erase(m_vertices.begin() + vertex2);
    m_vertices.erase(m_vertices.begin() + vertex3);
}
void ConvexHull::AddVertex(const Vector3& newVertex)
{
    // Add vertex to hull:
    m_hull.push_back(newVertex);
    // Remove all triangle it can see, construct polygonEdge of hole, and remove bad vertices:
    std::vector<Vector3Int> badTriangles = RemoveBadTriangles(newVertex);
    std::vector<Vector2Int> polygonEdges = PolygonEdges(badTriangles);
    RemoveBadVertices(badTriangles, polygonEdges);
    // connect new vertex to polygonEdge:
    center = Vector3::GetCenter(m_hull);
    for(const Vector2Int& edge : polygonEdges)
    {
        Vector3Int newTriangle = OrientedTriangle(m_hull.size() - 1, edge[0], edge[1]);
        m_triangles.push_back(newTriangle);
    }
}
bool ConvexHull::CanSeeTriangle(Vector3 vertex, Vector3Int triangle)
{
    Vector3 v0 = m_hull[triangle[0]];
    Vector3 v1 = m_hull[triangle[1]];
    Vector3 v2 = m_hull[triangle[2]];
    Vector3 n = Vector3::Cross(v1 - v0, v2 - v0).Normalized();
    Vector3 triangleCenter = (v0 + v1 + v2) / 3.0;
    return Vector3::Dot(n, triangleCenter - vertex) < 0;
}
Vector3Int ConvexHull::OrientedTriangle(int a, int b, int c)
{
    Vector3Int triangle(a,b,c);
    if (!CanSeeTriangle(center, triangle))
        return triangle;
    else
        return Vector3Int(a, c, b);
}
bool ConvexHull::InsideHull(Vector3 vertex)
{
    for(Vector3Int& triangle : m_triangles)
    {
        if (CanSeeTriangle(vertex, triangle))
            return false;
    }
    return true;
}

std::vector<Vector3Int> ConvexHull::RemoveBadTriangles(Vector3 vertex)
{
    std::vector<Vector3Int> badTriangles;
    for(int i=m_triangles.size()-1; i>-1; i--)
    {
        Vector3Int triangle = m_triangles[i];
        if (CanSeeTriangle(vertex, triangle))
        {
            badTriangles.push_back(triangle);
            m_triangles.erase(m_triangles.begin() + i);
        }
    }
    return badTriangles;
}
std::vector<Vector2Int> ConvexHull::PolygonEdges(const std::vector<Vector3Int>& badTriangles)
{
    // Get all edges, including duplicates:
    std::vector<Vector2Int> allEdges;
    for(const Vector3Int& triangle : badTriangles)
    {
        allEdges.push_back(Vector2Int(triangle[0], triangle[1]));
        allEdges.push_back(Vector2Int(triangle[1], triangle[2]));
        allEdges.push_back(Vector2Int(triangle[2], triangle[0]));
    }
    // Only pick edges that are unique (including edge inversion).
    std::vector<Vector2Int> polygonEdges;
    for(int i=0; i<allEdges.size(); i++)
    {
        Vector2Int currentEdge = allEdges[i];
        bool isUnique = true;
        for(int j=0; j<allEdges.size(); j++)
        {
            Vector2Int otherEdge = allEdges[j];
            if(i == j)
                continue;
            if((currentEdge[0] == otherEdge[0] && currentEdge[1] == otherEdge[1])
            || (currentEdge[0] == otherEdge[1] && currentEdge[1] == otherEdge[0]))
            {
                isUnique = false;
                break;
            }
        }
        if(isUnique)
            polygonEdges.push_back(currentEdge);
    }
    return polygonEdges;
}
void ConvexHull::RemoveBadVertices(const std::vector<Vector3Int>& badTriangles, const std::vector<Vector2Int>& polygonEdges)
{
    // use set to keep badVertices unique.
    std::set<int> badVertices;
    std::vector<int> badTrianglesUnrolled;
    for(const Vector3Int& triangle : badTriangles)
    {
        badTrianglesUnrolled.push_back(triangle[0]);
        badTrianglesUnrolled.push_back(triangle[1]);
        badTrianglesUnrolled.push_back(triangle[2]);
    }
    std::vector<int> polygonEdgesUnrolled;
    for(const Vector2Int& edge : polygonEdges)
    {
        polygonEdgesUnrolled.push_back(edge[0]);
        polygonEdgesUnrolled.push_back(edge[1]);
    }
    for(int triangleIndex : badTrianglesUnrolled)
    {
        bool goodIndex = false;
        for(int polygonIndex : polygonEdgesUnrolled)
            if (triangleIndex == polygonIndex)
            {
                goodIndex = true;
                break;
            }
        if (!goodIndex)
            badVertices.insert(triangleIndex);
    }
    // Remove badVertices in reversed order:
    for (auto rit = badVertices.rbegin(); rit != badVertices.rend(); rit++)
        m_vertices.erase(m_vertices.begin() + *rit);
}