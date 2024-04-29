#include "InterpolationGrid.h"

// Constructor:
InterpolationGrid::InterpolationGrid(size_t nTh, size_t nPh, const Stencil &stencil)
    : nTh(nTh), nPh(nPh), nDir(nTh * nPh)
{
    // Precalculate barycentric neighbours and weights on interpolation grid:
    barycentricNeighbours.reserve(nDir);
    barycentricWeights.reserve(nDir);
    for (size_t j = 0; j < nPh; j++)
        for (size_t i = 0; i < nTh; i++)
        {
            Vector3 p = Cv3(i, j);
            std::tuple<Vector3Int, Vector3> neighboursAndWeights = stencil.BarycentricNeighboursAndWeights(p);
            barycentricNeighbours.push_back(std::get<0>(neighboursAndWeights));
            barycentricWeights.push_back(std::get<1>(neighboursAndWeights));
        }

    // Precalculate voronoi neighbours and weights on interpolation grid:
    voronoiNeighbours.reserve(10 * nDir);
    voronoiWeights.reserve(10 * nDir);
    for (size_t j = 0; j < nPh; j++)
        for (size_t i = 0; i < nTh; i++)
        {
            Vector3 p = Cv3(i, j);
            std::tuple<std::vector<size_t>, std::vector<double>> neighboursAndWeights = stencil.VoronoiNeighboursAndWeights(p);
            voronoiNeighbours.AddRow(std::get<0>(neighboursAndWeights));
            voronoiWeights.AddRow(std::get<1>(neighboursAndWeights));
        }

    voronoiNeighbours.shrink_to_fit();
    voronoiWeights.shrink_to_fit();
}

// Grid operations:
double InterpolationGrid::Theta(size_t i) const
{
    return M_PI * i / (nTh - 1.0);
}
double InterpolationGrid::Phi(size_t j) const
{
    return 2.0 * M_PI * j / (double)nPh;
}
double InterpolationGrid::i(double theta) const
{
    return std::clamp((nTh - 1.0) * theta / M_PI, 1e-8, nTh - 1.0 - 1e-8);
}
double InterpolationGrid::j(double phi) const
{
    return std::clamp(nPh * phi / (2.0 * M_PI), 1e-8, nPh - 1e-8);
}
double InterpolationGrid::i(const Tensor3 &c) const
{
    return i(c.Theta());
}
double InterpolationGrid::j(const Tensor3 &c) const
{
    return std::fmod((j(c.Phi()) + nPh), nPh);
}
size_t InterpolationGrid::Index(size_t i, size_t j) const
{
    return i + j * nTh;
}
Tensor3 InterpolationGrid::Ct3(size_t i, size_t j) const
{
    double theta = Theta(i);
    double phi = Phi(j);
    return Tensor3(MySin(theta) * MyCos(phi), MySin(theta) * MySin(phi), MyCos(theta));
}
Vector3 InterpolationGrid::Cv3(size_t i, size_t j) const
{
    double theta = Theta(i);
    double phi = Phi(j);
    return Vector3(MySin(theta) * MyCos(phi), MySin(theta) * MySin(phi), MyCos(theta));
}

// Getters:
Vector3Int InterpolationGrid::BarycentricNeighbours(size_t i, size_t j) const
{
    int ij = Index(i, j);
    return barycentricNeighbours[ij];
}
Vector3 InterpolationGrid::BarycentricWeights(size_t i, size_t j) const
{
    int ij = Index(i, j);
    return barycentricWeights[ij];
}
std::span<const size_t> InterpolationGrid::VoronoiNeighbours(size_t i, size_t j) const
{
    int ij = Index(i, j);
    return voronoiNeighbours[ij];
}
std::span<const double> InterpolationGrid::VoronoiWeights(size_t i, size_t j) const
{
    int ij = Index(i, j);
    return voronoiWeights[ij];
}