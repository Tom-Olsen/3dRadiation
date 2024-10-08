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
    if (PointToPointDistance(Cv3(cellIndex), p) < 1e-8)
    {
        std::span<const size_t> naturalNeighboursSpan = voronoiNeighbours[cellIndex];
        naturalNeighbours = std::vector<size_t>(naturalNeighboursSpan.begin(), naturalNeighboursSpan.end());
        pointsInsideCell.clear();

        std::span<const Vector3> voronoiCellSpan = VoronoiCell(cellIndex);
        return std::vector<Vector3>(voronoiCellSpan.begin(), voronoiCellSpan.end());
    }

    // Declarations:
    std::vector<Vector3> virtualCell;
    Vector3 intersect;
    Vector3 rayOrigin;
    Vector3 prevPlaneNormal;
    Vector3 nextPrevPlaneNormal;
    for (size_t i = 0; i < 100; i++) // should be at most 6 itterations!
    {
        // Get current cell:
        std::span<const Vector3> cell = VoronoiCell(cellIndex);
        Vector3 cellCenter = Cv3(cellIndex);

        // Prepare line for LinePlaneIntersection:
        Vector3 dir = (cellCenter - p) / 2.0;
        Vector3 mid = (p + dir).Normalized();
        Vector3 perp = Vector3::Cross(mid, dir);
        if (i == 0)
            rayOrigin = mid;

        // Find intersection with the edge of the cell that is hit by perp head on:
        size_t index;
        double tMin = 1e10;
        for (size_t j = 0; j < cell.size(); j++)
        {
            Vector3 cellWallDir = cell[j] - cell[(j + 1) % cell.size()];
            Vector3 planeNormal = Vector3::Cross(cell[j], cellWallDir);

            // Skip origin edge (j=0 this will alway be false):
            if (IsParallel(prevPlaneNormal, planeNormal))
                continue;

            double t;
            Vector3 intersectionTemp = LinePlaneIntersection(rayOrigin, perp, cell[j], planeNormal, t);
            double dist = PointToPointDistance(rayOrigin, intersectionTemp);
            if (1e-8 < t && t < tMin)
            {
                tMin = t;
                index = j;
                intersect = intersectionTemp;
                nextPrevPlaneNormal = planeNormal;
            }
        }
        // Add to virtual cell:
        prevPlaneNormal = nextPrevPlaneNormal;
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
bool ContainsSimilar(const std::vector<Vector3> &list, Vector3 &p, double epsilon = 1e-8)
{
    for (const Vector3 &q : list)
        if ((p - q).Norm() < epsilon)
            return true;
    return false;
}
void Stencil::AddGhostDirections()
{
    double zThreshold0 = MyCos(refinement0Threshold * M_PI);
    double zThreshold1 = MyCos(refinement1Threshold * M_PI);
    double zThreshold2 = MyCos(refinement2Threshold * M_PI);
    std::vector<Vector3Int> triangles = mesh.GetTriangles();
    std::vector<Vector3> ghostDirections;

    // vertices:
    // for(int d=0; d<nDir; d++)
        // std::cout << "(" << Format(Cx(d),6) << "," << Format(Cy(d),6) << "," << Format(Cz(d),6) << "," << Format(Theta(d),6) << "," << Format(Phi(d),6) << ")" << std::endl;

    // Calculate extra refinement ghost directions:
    for (int t = 0; t < triangles.size(); t++)
    {
        Vector3Int triangle = triangles[t];
        Vector3 a = Cv3(triangle[0]);
        Vector3 b = Cv3(triangle[1]);
        Vector3 c = Cv3(triangle[2]);
        Vector3 p = (a + b + c).Normalized();
        double z = p[2];
        
        // triangles:
        // std::cout << "(" << triangle[0] << "," << triangle[1] << "," << triangle[2] << ")" << std::endl;
        // centers:
        // std::cout << "(" << Format(p[0],6) << "," << Format(p[1],6) << "," << Format(p[2],6) << "," << Format(p.Theta(),6) << "," << Format(p.Phi(),6) << ")" << std::endl;

        if (z > zThreshold2)
        {
            // focal centers:
            // std::cout << "(" << Format(p[0],6) << "," << Format(p[1],6) << "," << Format(p[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(p[0],6) << "," << Format(p[1],6) << "," << Format(p[2],6) << "," << Format(p.Theta(),6) << "," << Format(p.Phi(),6) << ")" << std::endl;
            
            Vector3 p0 = (a + b).Normalized();
            Vector3 p1 = (b + c).Normalized();
            Vector3 p2 = (c + a).Normalized();
            Vector3 q0 = (p0 + p1 + b).Normalized();
            Vector3 q1 = (p1 + p2 + c).Normalized();
            Vector3 q2 = (p2 + p0 + a).Normalized();

            // ghost directions:
            // std::cout << "(" << Format( p[0],6) << "," << Format( p[1],6) << "," << Format( p[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(p0[0],6) << "," << Format(p0[1],6) << "," << Format(p0[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(p1[0],6) << "," << Format(p1[1],6) << "," << Format(p1[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(p2[0],6) << "," << Format(p2[1],6) << "," << Format(p2[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(q0[0],6) << "," << Format(q0[1],6) << "," << Format(q0[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(q1[0],6) << "," << Format(q1[1],6) << "," << Format(q1[2],6) << ")" << std::endl;
            // std::cout << "(" << Format(q2[0],6) << "," << Format(q2[1],6) << "," << Format(q2[2],6) << ")" << std::endl;

            if (!ContainsSimilar(ghostDirections, p,  1e-4))
                ghostDirections.push_back(p);
            if (!ContainsSimilar(ghostDirections, p0, 1e-4))
                ghostDirections.push_back(p0);
            if (!ContainsSimilar(ghostDirections, p1, 1e-4))
                ghostDirections.push_back(p1);
            if (!ContainsSimilar(ghostDirections, p2, 1e-4))
                ghostDirections.push_back(p2);
            if (!ContainsSimilar(ghostDirections, q0, 1e-4))
                ghostDirections.push_back(q0);
            if (!ContainsSimilar(ghostDirections, q1, 1e-4))
                ghostDirections.push_back(q1);
            if (!ContainsSimilar(ghostDirections, q2, 1e-4))
                ghostDirections.push_back(q2);
        }
        else if (z > zThreshold1)
        {
            Vector3 p0 = (a + b).Normalized();
            Vector3 p1 = (b + c).Normalized();
            Vector3 p2 = (c + a).Normalized();
            if (!ContainsSimilar(ghostDirections, p0, 1e-4))
                ghostDirections.push_back(p0);
            if (!ContainsSimilar(ghostDirections, p1, 1e-4))
                ghostDirections.push_back(p1);
            if (!ContainsSimilar(ghostDirections, p2, 1e-4))
                ghostDirections.push_back(p2);
        }
        else if (z > zThreshold0)
        {
            if (!ContainsSimilar(ghostDirections, p, 1e-4))
                ghostDirections.push_back(p);
        }
    }

    // Add ghost directions to velocity stencil:
    for (const Vector3 &p : ghostDirections)
    {
        w.push_back(0);
        cx.push_back(p[0]);
        cy.push_back(p[1]);
        cz.push_back(p[2]);
        theta.push_back(p.Theta());
        phi.push_back(p.Phi());
    }
    nGhost = ghostDirections.size();
    nDir += nGhost;
}
void Stencil::SortDirections()
{
    // Combine arrays into a vector of tuples
    std::vector<std::tuple<double, double, double, double, double, double>> combinedArrays;
    for (int i = 0; i < nDir; ++i)
        combinedArrays.emplace_back(w[i], theta[i], phi[i], cx[i], cy[i], cz[i]);

    // Sort the vector of tuples based on theta and then phi
    std::sort(combinedArrays.begin(), combinedArrays.end(), 
        [](const auto &a, const auto &b)
        {
            // Round theta and phi values
            double roundedThetaA  = Round(std::get<1>(a), 6);
            double roundedPhiA    = Round(std::get<2>(a), 6);
            double roundedThetaB  = Round(std::get<1>(b), 6);
            double roundedPhiB    = Round(std::get<2>(b), 6);

            // Compare rounded theta and phi
            if (roundedThetaA != roundedThetaB)
                return roundedThetaA < roundedThetaB;
            else
                return roundedPhiA < roundedPhiB;
        }
    );

    // Write sorted values back into original arrays
    for (int i = 0; i < nDir; ++i)
    {
        w[i] = std::get<0>(combinedArrays[i]);
        theta[i] = std::get<1>(combinedArrays[i]);
        phi[i] = std::get<2>(combinedArrays[i]);
        cx[i] = std::get<3>(combinedArrays[i]);
        cy[i] = std::get<4>(combinedArrays[i]);
        cz[i] = std::get<5>(combinedArrays[i]);
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
        voronoiCell.reserve(triangles.size());
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
        neighbours.shrink_to_fit();
        voronoiNeighbours.AddRow(neighbours);
    }
    voronoiNeighbours.shrink_to_fit();
}

void Stencil::PopulateLookUpTable()
{
    sigmaMax = MaxSigma();
    double deltaSigma = 0.1;
    for(double sigma = 0.0; sigma <= sigmaMax; sigma += deltaSigma)
    {
        double flux = (sigma * cosh(sigma) - sinh(sigma)) / (sigma * sinh(sigma));
        if (sigma == 0.0)
            sigmaOfRelativeFlux.Add(0, 0);
        else
            sigmaOfRelativeFlux.Add(flux, sigma);

        if (sigmaMax - sigma < deltaSigma)
            deltaSigma = sigmaMax - sigma;
        if (sigmaMax - sigma < 1e-8)
            break;
    }
    relativeFluxMax = sigmaOfRelativeFlux.xValues[sigmaOfRelativeFlux.size - 1];
}
double Stencil::MaxSigma()
{
    RealBuffer I;
    I.resize(nDir);
    int nThTestGrid = 100;
    int nPhTestGrid = 200;
    double sigma;
    double deltaSigma = 1;

    for(sigma = 1; sigma < 500; sigma += deltaSigma)
    {
        for (int d = 0; d < nDir; d++)
            I[d] = Intensity(sigma, 1.0, Theta(d));

        double averageError = 0;
        for (int d1 = 0; d1 < nPhTestGrid; d1++)
            for (int d0 = 0; d0 < nThTestGrid; d0++)
            {
                int d = d0 + d1 * nThTestGrid;
                // Only look at theta=0->0.1 M_PI disk around north pole:
                double theta = 0.1 * M_PI * d0 / (nThTestGrid - 1.0);
                double phi = 2.0 * M_PI * d1 / (double)nPhTestGrid;
                Vector3 dir = Vector3(MySin(theta) * MyCos(phi), MySin(theta) * MySin(phi), MyCos(theta));
    
                // Analytic intensity:
                double analyticValue = Intensity(sigma, 1.0, theta);
    
                // Voronoi interpolated intensity:
                std::tuple<std::vector<size_t>, std::vector<double>> neighboursAndWeights = VoronoiNeighboursAndWeights(dir);
                std::span<const size_t> neighbours = std::get<0>(neighboursAndWeights);
                std::span<const double> weights = std::get<1>(neighboursAndWeights);
                double interpolatetValue = 0;
                for (size_t p = 0; p < weights.size(); p++)
                    interpolatetValue += weights[p] * I[neighbours[p]];
    
                // Error:
                double error = std::abs((analyticValue - interpolatetValue) / analyticValue);
                averageError += error;
            }
        averageError /= (nThTestGrid * nPhTestGrid);

        // Check if error is too small:
        if (averageError >= maxInterpolationError)
        {
            // Roll back sigma and reduce deltaSigma:
            sigma -= deltaSigma;
            deltaSigma /= 10;
            
            // Stop after 6 sigits of accuracy:
            if (deltaSigma <= sigmaMaxSearchAccuracy)
                return sigma;
        }
    }

    ExitOnError("MaxSigma: Could not find sigma with error below " + std::to_string(maxInterpolationError));
    return 0;
}

// Debugging:
void Stencil::Print() const
{
    std::cout << name << ":\n";
    std::cout << "           nDir: " << nDir << "\n";
    std::cout << "         nGhost: " << nGhost << "\n";
    std::cout << "       sigmaMax: " << sigmaMax << "\n";
    std::cout << "relativeFluxMax: " << relativeFluxMax << "\n";
}
void Stencil::PrintAll() const
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
LebedevStencil::LebedevStencil(size_t nOrder, double refinement0Threshold, double refinement1Threshold, double refinement2Threshold)
{
    this->name = "Lebedev" + std::to_string(nOrder) + "_" + std::to_string(refinement0Threshold) + "_" + std::to_string(refinement1Threshold) + "_" + FormatNoSignSpace(refinement2Threshold);
    this->nDir = 0;
    this->nOrder = nOrder;
    this->nCoefficients = IntegerPow<2>((nOrder + 1) / 2);
    this->refinement0Threshold = refinement0Threshold;
    this->refinement1Threshold = refinement1Threshold;
    this->refinement2Threshold = refinement2Threshold;

    std::ifstream inputFile("../stencils/LebedevStencil/LebedevStencil" + std::to_string(nOrder));
    if (!inputFile.is_open())
        ExitOnError("File '/stencils/LebedevStencil/LebedevStencil" + std::to_string(nOrder) + "' is missing.");

    // Read file:
    std::string line;
    int i = 0;
    while (std::getline(inputFile, line))
    {
        // Skip empty lines and comments:
        if (line.empty() || line[0] == '#')
            continue;
        else if (this->nDir == 0)
        {
            // Read count:
            this->nDir = std::stoi(line);
            AllocateBuffers();
        }
        else
        {
            // Parse data:
            std::istringstream iss(line);
            std::string lineStream;

            std::vector<double> row;
            while (getline(iss, lineStream, ','))
                row.push_back(stold(lineStream)); // convert to double

            w[i] = row[0];
            cx[i] = row[1];
            cy[i] = row[2];
            cz[i] = row[3];
            theta[i] = row[4];
            phi[i] = row[5];
            i++;
        }
    }
    inputFile.close();
    
    SortDirections();
    InitializeMesh();
    AddGhostDirections();
    SortDirections();
    InitializeMesh();
    InitializeConnectedTriangles();
    InitializeVoronoiCells();
    InitializeVoronoiNeighbours();
    PopulateLookUpTable();
}
// ---------------------------------------------------------------------

// ----------------------- GaussLegendreStencil ------------------------
GaussLegendreStencil::GaussLegendreStencil(size_t nOrder)
{
    this->name = "GaussLegendre" + std::to_string(nOrder);
    this->nDir = 0;
    this->nOrder = nOrder;
    this->nCoefficients = IntegerPow<2>((nOrder + 1) / 2);
    this->refinement0Threshold = 0;
    this->refinement1Threshold = 0;
    this->refinement2Threshold = 0;

    std::ifstream inputFile("../stencils/GaussLegendreStencil/GaussLegendreStencil" + std::to_string(nOrder));
    if (!inputFile.is_open())
        ExitOnError("File '/stencils/GaussLegendreStencil/GaussLegendreStencil" + std::to_string(nOrder) + "' is missing.");

    std::string line;
    int i = 0;
    while (std::getline(inputFile, line))
    {
        // Skip empty lines and comments:
        if (line.empty() || line[0] == '#')
            continue;
        else if (this->nDir == 0)
        {
            // Read count:
            this->nDir = std::stoi(line);
            AllocateBuffers();
        }
        else
        {
            // Parse data:
            std::istringstream iss(line);
            std::string lineStream;

            std::vector<double> row;
            while (getline(iss, lineStream, ','))
                row.push_back(stod(lineStream)); // convert to double

            w[i] = row[0];
            cx[i] = row[1];
            cy[i] = row[2];
            cz[i] = row[3];
            theta[i] = row[4];
            phi[i] = row[5];
            i++;
        }
    }
    inputFile.close();

    SortDirections();
    InitializeMesh();
    InitializeConnectedTriangles();
    InitializeVoronoiCells();
    InitializeVoronoiNeighbours();
    PopulateLookUpTable();
}
// ---------------------------------------------------------------------