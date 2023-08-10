#include "Stencil.h"

// ------------------------------- Stencil -------------------------------
// Getters:
double Stencil::W(size_t d) const
{
    return w[d];
}
double Stencil::Cx(size_t d) const
{
    return cx[d];
}
double Stencil::Cy(size_t d) const
{
    return cy[d];
}
double Stencil::Cz(size_t d) const
{
    return cz[d];
}
double Stencil::Theta(size_t d) const
{
    return theta[d];
}
double Stencil::Phi(size_t d) const
{
    return phi[d];
}
Tensor3 Stencil::Ct3(size_t d) const
{
    return Tensor3(Cx(d), Cy(d), Cz(d));
}
Vector3 Stencil::Cv3(size_t d) const
{
    return Vector3(Cx(d), Cy(d), Cz(d));
}
std::span<const Vector3Int> Stencil::ConnectedTriangles(size_t d) const
{
    return connectedTriangles[d];
}
std::span<const Vector3> Stencil::VoronoiCell(size_t d) const
{
    return voronoiCells[d];
}
std::span<const size_t> Stencil::VoronoiNeighbours(size_t d) const
{
    return voronoiNeighbours[d];
}

// Expensive Getters:
// Get voronoi neighbours and weights of an abitrary point p.
std::tuple<std::vector<size_t>, std::vector<double>> Stencil::VoronoiNeighboursAndWeights(const Vector3 &p) const
{
    std::vector<size_t> naturalNeighbours;
    std::vector<Vector3> pointsInsideCell;
    std::vector<Vector3> virtualCell = VirtualVoronoiCellOf(p, naturalNeighbours, pointsInsideCell);
    if (pointsInsideCell.size() == 0)
    {
        size_t index = NearestNeighbour(p);
        std::tuple<std::vector<size_t>, std::vector<double>> output;
        std::get<0>(output) = std::vector<size_t>{index};
        std::get<1>(output) = std::vector<double>{1.0};
        return output;
    }
    std::vector<double> weights = VoronoiWeights(virtualCell, naturalNeighbours, pointsInsideCell);
    return {naturalNeighbours, weights};
}
// Get barycentric neighbours and weights of an abitrary point p.
std::tuple<Vector3Int, Vector3> Stencil::BarycentricNeighboursAndWeights(const Vector3 &p) const
{
    Vector3 weights;
    Vector3 rayOrigin(0, 0, 0);
    std::vector<Vector3> vertices = mesh.GetVertices();
    std::vector<Vector3Int> triangles = mesh.GetTriangles();
    for (Vector3Int triangle : triangles)
    {
        Vector3 v0 = vertices[triangle[0]];
        Vector3 v1 = vertices[triangle[1]];
        Vector3 v2 = vertices[triangle[2]];
        if (BarycentricWeights(rayOrigin, p, v0, v1, v2, weights))
            return {triangle, weights};
    }
    return {Vector3Int(), Vector3()};
}

// Internal helper methods:
// Index of nearest stencil.C(?) for arbitrary point p.
size_t Stencil::NearestNeighbour(const Vector3 &p) const
{
    size_t index = -1;
    double dotMax = -1;
    for (size_t d = 0; d < nDir; d++)
    {
        double dot = Vector3::Dot(p, Cv3(d));
        if (dot > dotMax)
        {
            dotMax = dot;
            index = d;
        }
    }
    return index;
}
// Creates virtual voronoi cell of arbitrary point p. This includes the natural neighbours and voronoi points that are inside the virual cell.
std::vector<Vector3> Stencil::VirtualVoronoiCellOf(const Vector3 &p, std::vector<size_t> &naturalNeighbours, std::vector<Vector3> &pointsInsideCell) const
{
    // Initial cell given by nearest direction:
    size_t startCellIndex = NearestNeighbour(p);
    size_t cellIndex = startCellIndex;

    // Check if p overlaps with central vertex of cell:
    if (Vector3::Dot(Cv3(cellIndex) - p, Cv3(cellIndex) - p) < 1e-8)
    {
        std::span<const Vector3> temp = VoronoiCell(cellIndex);
        return std::vector<Vector3>(temp.begin(), temp.end());
    }

    // Declarations:
    std::vector<Vector3> virtualCell;
    Vector3 intersect;
    Vector3 intersectTemp;
    Vector3 rayOrigin;
    for (size_t i = 0; i < 100; i++) // should be at most 6 itterations!
    {
        // Get current localized cell:
        std::span<const Vector3> cell = VoronoiCell(cellIndex);
        Vector3 cellCenter = Cv3(cellIndex);

        // Prepare geometry for ray-lineSegment intersection:
        Vector3 dir = (cellCenter - p) / 2.0;
        Vector3 mid = (p + dir).Normalized();
        Vector3 perp = Vector3::Cross(mid, dir);
        if (i == 0)
            rayOrigin = mid;

        // Find intersection with the edge of the cell that is hit by perp head on:
        size_t index;
        double tMin = 1e10;
        double t;
        for (size_t j = 0; j < cell.size(); j++)
        {
            Vector3 cellWallDir = cell[j] - cell[(j + 1) % cell.size()];
            Vector3 planeNormal = Vector3::Cross(cell[j], cellWallDir);

            if (RayPlaneIntersection(rayOrigin, perp, cell[j], planeNormal, intersectTemp, t))
            {
                if (t < tMin)
                {
                    tMin = t;
                    index = j;
                    intersect = intersectTemp;
                }
            }
        }
        // Add to virtual cell:
        intersect = intersect.Normalized();
        virtualCell.push_back(intersect);

        // Add vertex which will end up inside the virtual cell:
        if (!Contains(pointsInsideCell, cell[index]))
            pointsInsideCell.push_back(cell[index]);

        // Prepare next cell: (only works because of ordering of cells and orientation of virtual cell)
        rayOrigin = intersect;
        cellIndex = VoronoiNeighbours(cellIndex)[index];
        naturalNeighbours.push_back(cellIndex);

        // End loop once we are back at initial cell:
        if (cellIndex == startCellIndex)
            break;
    }
    return virtualCell;
}
// Get voronoi weights of virtual voronoi cell.
std::vector<double> Stencil::VoronoiWeights(const std::vector<Vector3> &virtualCell, const std::vector<size_t> &naturalNeighbours, const std::vector<Vector3> &pointsInsideCell) const
{
    double totalArea = 0.0;
    std::vector<double> weights;
    weights.reserve(naturalNeighbours.size());
    for (size_t i = 0; i < naturalNeighbours.size(); i++)
    {
        // Construct intersection polygon of virtualCell and cell:
        std::span<const Vector3> cell = VoronoiCell(naturalNeighbours[i]);
        std::vector<Vector3> points;
        points.push_back(virtualCell[i]);
        points.push_back(virtualCell[(i + 1) % naturalNeighbours.size()]);
        for (Vector3 p : pointsInsideCell)
            if (Contains(cell, p))
                points.push_back(p);

        // Get center of polygon:
        Vector3 center = Vector3(0, 0, 0);
        for (Vector3 p : points)
            center += p;
        center = center.Normalized();

        // Sort polygon around center:
        points = SortPointsOnSphereAroundCenter(points);

        // Slice polycon into triangles and determine area:
        double area = 0.0;
        for (size_t j = 0; j < points.size(); j++)
            area += SphericalTriangleArea(points[j], points[(j + 1) % points.size()], center);
        weights.push_back(area);
        totalArea += area;
    }
    // Normalize weights:
    for (size_t i = 0; i < weights.size(); i++)
        weights[i] /= totalArea;
    return weights;
}

// Initialization:
void Stencil::AllocateBuffers()
{
    w.resize(nDir);
    cx.resize(nDir);
    cy.resize(nDir);
    cz.resize(nDir);
    theta.resize(nDir);
    phi.resize(nDir);
}
void Stencil::AddGhostDirections()
{
    size_t index = nDir - nGhost;
    for (size_t i = 0; i < nRings; i++)
    {
        double ringTheta = thetaGhost * (i + 0.5) / nRings;
        for (size_t j = 0; j < nRing0 * (i + 1); j++)
        {
            w[index] = 0.0;
            theta[index] = ringTheta;
            phi[index] = 2.0 * M_PI * j / (nRing0 * (i + 1.0));
            cx[index] = MySin(theta[index]) * MyCos(phi[index]);
            cy[index] = MySin(theta[index]) * MySin(phi[index]);
            cz[index] = MyCos(theta[index]);
            index++;
        }
    }
}
void Stencil::SortDirections()
{
    size_t index[nDir];
    for (size_t i = 0; i < nDir; i++)
        index[i] = i;

    // Sort index array based on custom compare function:
    std::sort(index, index + nDir, [this](size_t i, size_t j)
              {
        if(theta[i] != theta[j])
            return theta[i] < theta[j];
        else
            return phi[i] < phi[j]; });

    // Create temporary arrays to hold sorted values:
    double sorted_w[nDir];
    double sorted_cx[nDir];
    double sorted_cy[nDir];
    double sorted_cz[nDir];
    double sorted_theta[nDir];
    double sorted_phi[nDir];

    // Copy values from original arrays to temporary arrays:
    for (size_t i = 0; i < nDir; i++)
    {
        sorted_w[i] = w[index[i]];
        sorted_cx[i] = cx[index[i]];
        sorted_cy[i] = cy[index[i]];
        sorted_cz[i] = cz[index[i]];
        sorted_theta[i] = theta[index[i]];
        sorted_phi[i] = phi[index[i]];
    }

    // Overwrite original arrays with sorted values:
    for (size_t i = 0; i < nDir; i++)
    {
        w[i] = sorted_w[i];
        cx[i] = sorted_cx[i];
        cy[i] = sorted_cy[i];
        cz[i] = sorted_cz[i];
        theta[i] = sorted_theta[i];
        phi[i] = sorted_phi[i];
    }
}
void Stencil::InitializeMesh()
{
    // Setup vertices vector for convex hull triangulation:
    std::vector<Vector3> vertices;
    vertices.reserve(nDir);
    for (size_t d = 0; d < nDir; d++)
        vertices.push_back(Cv3(d));

    // Triangulate vertices:
    ConvexHull convexHull(vertices);
    convexHull.OriginalOrdering(vertices);
    mesh = convexHull.GetMesh();
}
void Stencil::InitializeConnectedTriangles()
{
    // At most 10 natural neighbours (very generous estimate):
    connectedTriangles.reserve(10 * nDir);
    std::vector<Vector3Int> allTriangles = mesh.GetTriangles();
    for (size_t i = 0; i < nDir; i++)
    {
        std::vector<Vector3Int> triangles;
        triangles.reserve(10);
        for (size_t j = 0; j < allTriangles.size(); j++)
        {
            Vector3Int triangle = allTriangles[j];
            if (triangle[0] == i || triangle[1] == i || triangle[2] == i)
                triangles.push_back(triangle);
        }
        // no triangles.shrink_to_fit() becuase after copying triangles will be discarded.
        connectedTriangles.AddRow(triangles);
    }
    connectedTriangles.shrink_to_fit();
}
void Stencil::InitializeVoronoiCells()
{
    // At most 10 natural neighbours (very generous estimate):
    voronoiCells.reserve(10 * nDir);
    for (size_t d = 0; d < nDir; d++)
    {
        // Get circum center of all triangles connected to a vertex:
        std::span<const Vector3Int> triangles = ConnectedTriangles(d);
        std::vector<Vector3> voronoiCell;
        voronoiCell.reserve(triangles.size()); // no need for shirnk_to_fit.
        for (const Vector3Int &triangle : triangles)
        {
            Vector3 a = Cv3(triangle[0]);
            Vector3 b = Cv3(triangle[1]);
            Vector3 c = Cv3(triangle[2]);
            Vector3 center = CircumCenter(a, b, c).Normalized();
            if (!Contains(voronoiCell, center))
                voronoiCell.push_back(center);
        }
        voronoiCell.shrink_to_fit();
        voronoiCell = SortPointsOnSphereAroundCenter(voronoiCell);
        voronoiCells.AddRow(voronoiCell);
    }
    voronoiCells.shrink_to_fit();
}
void Stencil::InitializeVoronoiNeighbours()
{
    voronoiNeighbours.reserve(10 * nDir);
    for (size_t d = 0; d < nDir; d++)
    {
        std::span<const Vector3> cell = VoronoiCell(d);
        std::vector<size_t> neighbours;
        neighbours.reserve(cell.size());
        for (size_t i = 0; i < cell.size(); i++)
        {
            Vector3 p0 = cell[i];
            Vector3 p1 = cell[(i + 1) % cell.size()];
            for (size_t k = 0; k < nDir; k++)
            {
                if (d == k)
                    continue;
                std::span<const Vector3> otherCell = VoronoiCell(k);
                if (Contains(otherCell, p0) && Contains(otherCell, p1))
                    neighbours.push_back(k);
            }
        }
        voronoiNeighbours.AddRow(neighbours);
    }
    voronoiNeighbours.shrink_to_fit();
}

// Move to InterpolationGrid?
// void Stencil::InitializeNearestNeighbourOnGrid()
//{
//    interpolationGrid = InterpolationGrid(interpolationGridRes, 2 * interpolationGridRes, *this);
//    PARALLEL_FOR(2)
//    for(size_t j=0; j<interpolationGrid.nPh; j++)
//    for(size_t i=0; i<interpolationGrid.nTh; i++)
//    {
//        Tensor3 c = interpolationGrid.Ct3(i,j);
//        double dot0 = -1;
//        double dot1 = -1;
//        double dot2 = -1;
//        size_t index0 = -1;
//        size_t index1 = -1;
//        size_t index2 = -1;
//        for(size_t d=0; d<nDir; d++)
//        {
//            double dot = Tensor3::Dot(c,Ct3(d));
//            if(dot > dot0)
//            {
//                dot2 = dot1;
//                dot1 = dot0;
//                dot0 = dot;
//                index2 = index1;
//                index1 = index0;
//                index0 = d;
//            }
//            else if(dot > dot1)
//            {
//                dot2 = dot1;
//                dot1 = dot;
//                index2 = index1;
//                index1 = d;
//            }
//            else if(dot > dot2)
//            {
//                dot2 = dot;
//                index2 = d;
//            }
//        }
//        neighbour0OnGrid[interpolationGrid.Index(i,j)] = index0;
//        neighbour1OnGrid[interpolationGrid.Index(i,j)] = index1;
//        neighbour2OnGrid[interpolationGrid.Index(i,j)] = index2;
//    }
//}
// void Stencil::InitializeVoronoiInterpolationOnGrid()
//{
//    for(size_t j=0; j<interpolationGrid.nPh; j++)
//    for(size_t i=0; i<interpolationGrid.nTh; i++)
//    {
//        Tensor3 c = interpolationGrid.Ct3(i,j);
//        std::vector<size_t> neighbours;
//        std::vector<double> weights = VoronoiWeights(c,neighbours);
//        voronoiNeighboursOnGrid.AddRow(neighbours);
//        voronoiWeightsOnGrid.AddRow(weights);
//    }
//}

// Debugging:
void Stencil::Print() const
{
    std::cout << name << ":\n";
    std::cout << "        d,\t        w,\t    theta,\t      phi,\t       cx,\t       cy,\t       cz\n";
    for (size_t d = 0; d < nDir; d++)
    {
        std::cout << Format(d) << ",\t" << Format(W(d)) << ",\t" << Format(Theta(d)) << ",\t" << Format(Phi(d)) << ",\t";
        std::cout << Format(Cx(d)) << ",\t" << Format(Cy(d)) << ",\t" << Format(Cz(d)) << "\n";
    }
    std::cout << std::endl;
}
// -----------------------------------------------------------------------

// -------------------------- LebedevStencil ---------------------------
LebedevStencil::LebedevStencil(size_t nOrder, size_t nRings, size_t nRing0, double thetaGhost)
{
    this->name = "Lebedev" + std::to_string(nOrder) + "." + std::to_string(nRings) + "." + std::to_string(nRing0) + "_" + FormatNoSignSpace(thetaGhost / M_PI, 3) + "pi";
    this->nOrder = nOrder;
    this->nCoefficients = IntegerPow<2>((nOrder + 1) / 2);
    this->nRings = nRings;
    this->nRing0 = nRing0;
    this->nGhost = nRing0 * nRings * (nRings + 1) / 2;
    this->thetaGhost = thetaGhost;

    switch (nOrder)
    {
    case 3:
#include "../stencils/LebedevStencil/LebedevStencil3"
        break;
    case 5:
#include "../stencils/LebedevStencil/LebedevStencil5"
        break;
    case 7:
#include "../stencils/LebedevStencil/LebedevStencil7"
        break;
    case 9:
#include "../stencils/LebedevStencil/LebedevStencil9"
        break;
    case 11:
#include "../stencils/LebedevStencil/LebedevStencil11"
        break;
        //    case 13:
        // #include "../stencils/LebedevStencil/LebedevStencil13"
        //        break;
    case 15:
#include "../stencils/LebedevStencil/LebedevStencil15"
        break;
    case 17:
#include "../stencils/LebedevStencil/LebedevStencil17"
        break;
    case 19:
#include "../stencils/LebedevStencil/LebedevStencil19"
        break;
    case 21:
#include "../stencils/LebedevStencil/LebedevStencil21"
        break;
    case 23:
#include "../stencils/LebedevStencil/LebedevStencil23"
        break;
        //    case 25:
        // #include "../stencils/LebedevStencil/LebedevStencil25"
        //        break;
        //    case 27:
        // #include "../stencils/LebedevStencil/LebedevStencil27"
        //        break;
    case 29:
#include "../stencils/LebedevStencil/LebedevStencil29"
        break;
    case 31:
#include "../stencils/LebedevStencil/LebedevStencil31"
        break;
    case 35:
#include "../stencils/LebedevStencil/LebedevStencil35"
        break;
    case 41:
#include "../stencils/LebedevStencil/LebedevStencil41"
        break;
    case 47:
#include "../stencils/LebedevStencil/LebedevStencil47"
        break;
    case 53:
#include "../stencils/LebedevStencil/LebedevStencil53"
        break;
    case 59:
#include "../stencils/LebedevStencil/LebedevStencil59"
        break;
    default:
        ExitOnError("LebedevStencil nOrder=" + std::to_string(nOrder) + " invalid. See stencil/LebedevStencil for valid orders.");
        break;
    }

    AddGhostDirections();
    SortDirections();
    InitializeMesh();
    InitializeConnectedTriangles();
    InitializeVoronoiCells();
    InitializeVoronoiNeighbours();
    // InitializeNearestNeighbourOnGrid();
    // InitializeVoronoiInterpolationOnGrid();
}
// ---------------------------------------------------------------------

// ----------------------- GaussLegendreStencil ------------------------
GaussLegendreStencil::GaussLegendreStencil(size_t nOrder)
{
    this->name = "GaussLegendre" + std::to_string(nOrder);
    this->nOrder = nOrder;
    this->nCoefficients = IntegerPow<2>((nOrder + 1) / 2);
    this->nRings = 0;
    this->nRing0 = 0;
    this->nGhost = 0;
    this->thetaGhost = 0;

    switch (nOrder)
    {
    case 3:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil3"
        break;
    case 5:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil5"
        break;
    case 7:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil7"
        break;
    case 9:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil9"
        break;
    case 11:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil11"
        break;
    case 13:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil13"
        break;
    case 15:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil15"
        break;
    case 17:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil17"
        break;
    case 19:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil19"
        break;
    case 21:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil21"
        break;
    case 23:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil23"
        break;
    case 25:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil25"
        break;
    case 27:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil27"
        break;
    case 29:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil29"
        break;
    case 31:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil31"
        break;
    case 33:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil33"
        break;
    case 35:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil35"
        break;
    case 37:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil37"
        break;
    case 39:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil39"
        break;
    case 41:
#include "../stencils/GaussLegendreStencil/GaussLegendreStencil41"
        break;
    default:
        ExitOnError("GaussLegendreStencil nOrder=" + std::to_string(nOrder) + " invalid. See stencil/GaussLegendreStencil for valid orders.");
        break;
    }

    SortDirections();
    InitializeMesh();
    InitializeConnectedTriangles();
    InitializeVoronoiCells();
    InitializeVoronoiNeighbours();
    // InitializeNearestNeighbourOnGrid();
    // InitializeVoronoiInterpolationOnGrid();
}
// ---------------------------------------------------------------------