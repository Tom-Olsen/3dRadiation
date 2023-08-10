#ifndef __INCLUDE_GUARD_InterpolationGrid_h__
#define __INCLUDE_GUARD_InterpolationGrid_h__
#include "Stencil.h"

struct InterpolationGrid
{
public:
    // Public Properties:
    size_t nTh;
    size_t nPh;
    size_t nDir;

private:
    // Internal Buffers (data for each grid point):
    std::vector<Vector3Int> barycentricNeighbours;
    std::vector<Vector3> barycentricWeights;
    UnstructuredMatrix<size_t> voronoiNeighbours;
    UnstructuredMatrix<double> voronoiWeights;

public:
    // Constructor:
    InterpolationGrid() = default;
    InterpolationGrid(size_t nTh, size_t nPh, const Stencil &stencil);

    // Grid operations:
    double Theta(size_t i) const;
    double Phi(size_t j) const;
    double i(double theta) const;
    double j(double phi) const;
    double i(const Tensor3 &c) const;
    double j(const Tensor3 &c) const;
    size_t Index(size_t i, size_t j) const;
    Tensor3 Ct3(size_t i, size_t j) const;
    Vector3 Cv3(size_t i, size_t j) const;

    // Getters:
    Vector3Int BarycentricNeighbours(size_t i, size_t j) const;
    Vector3 BarycentricWeights(size_t i, size_t j) const;
    std::span<const size_t> VoronoiNeighbours(size_t i, size_t j) const;
    std::span<const double> VoronoiWeights(size_t i, size_t j) const;

public:
    // Debugging:
    void Print() const;
};
#endif //__INCLUDE_GUARD_InterpolationGrid_h__