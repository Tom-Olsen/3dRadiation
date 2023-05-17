#ifndef __INCLUDE_GUARD_Stencil_h__
#define __INCLUDE_GUARD_Stencil_h__
#include <vector>
#include <set>
#include <span>
#include "Utility.hh"
#include "TensorTypes.hh"
#include "ConvexHull.h"
#include "Interpolation.h"
#include "SphereGrid.h"
#include "AdvancedMath.h"



// Spherical Coordinates convention:
// theta € [0,pi]
//   phi € [0,2pi]
// theta =  0 => North Pole
// theta = pi => South Pole
// (th,phi) = (pi/2,    0) => (x,y,z) = ( 1, 0, 0)
// (th,phi) = (pi/2, pi/2) => (x,y,z) = ( 0, 1, 0)
// (th,phi) = (pi/2,   pi) => (x,y,z) = (-1, 0, 0)
// (th,phi) = (pi/2,3pi/2) => (x,y,z) = ( 0,-1, 0)
// (th,phi) = (pi/2,  2pi) => (x,y,z) = ( 1, 0, 0)



// Base Stencil class. This functionality is supported by all stencils.
struct Stencil
{
public: // Public Properties:
    std::string name;       // Name of the stencel, e.g. "Lebedev9"
    size_t nDir;            // Number of directions in stencil.
    size_t nOrder;          // Quadrature integration.
    size_t nCoefficients;   // Number of exact Spherical Harmonics, counted with flat index i = l * (l + 1) + m.
    Mesh mesh;
    SphereGrid sphereGrid;
private:
    int sphereGridRes = 1000;

protected: // Internal Buffers on directions:
    RealBuffer w;
    RealBuffer cx;
    RealBuffer cy;
    RealBuffer cz;
    RealBuffer theta;
    RealBuffer phi;
    
    UnstructuredMatrix<Vector3Int> connectedTriangles;
    UnstructuredMatrix<Vector3> voronoiCells;
    UnstructuredMatrix<size_t> voronoiNeighbours;

protected: // Internal Buffers on sphereGrid:
    SizeTBuffer neighbour0OnGrid;
    SizeTBuffer neighbour1OnGrid;
    SizeTBuffer neighbour2OnGrid;
public: // Internal Buffers on sphereGrid:
    UnstructuredMatrix<size_t> voronoiNeighboursOnGrid;
    UnstructuredMatrix<double> voronoiWeightsOnGrid;

public: // Getters:
    double W(size_t d) const;
    double Cx(size_t d) const;
    double Cy(size_t d) const;
    double Cz(size_t d) const;
    double Theta(size_t d) const;
    double Phi(size_t d) const;
    Tensor3 Ct3(size_t d) const;
    Vector3 Cv3(size_t d) const;

    size_t NearestNeighbour(const Tensor3& p) const;
    std::span<const Vector3Int> TrianglesConnectedTo(size_t d) const;
    std::span<const Vector3> VoronoiCellOf(size_t d) const;
    std::span<const size_t> VoronoiNeighbourOf(size_t d) const;

    Vector3 BarycentricWeights(const Tensor3& p, Vector3Int& triangle) const;
    std::vector<double> VoronoiWeights(const Tensor3& p, std::vector<size_t>& naturalNeighbours) const;

protected:
    std::vector<Vector3> VirtualVoronoiCellOf
    (const Tensor3& p, std::vector<size_t>& naturalNeighbours, std::vector<Vector3>& pointsInsideCell) const;
    std::vector<double> VoronoiWeights
    (const std::vector<Vector3>& virtualCell, const std::vector<size_t>& naturalNeighbours, const std::vector<Vector3>& pointsInsideCell) const;

protected: // Initialization:
    void SetCoefficientCount();
    void AllocateBuffers();
    void SortDirections();
    void InitializeMesh();
    void InitializeConnectedTriangles();
    void InitializeVoronoiCells();
    void InitializeVoronoiNeighbours();
    void InitializeNearestNeighbourOnGrid();
    void InitializeVoronoiInterpolationOnGrid();

public: // Debugging:
    void Print() const;
};



// This class holds discretized velocities and the corresponding weights.
// The distribution of the velocity vectors is given by the spherical Lebedev quadarture.
// Exact Spherical Harmonic integration:
// LebedevStencil(1): Y00,
// LebedevStencil(3): Y00, Y1m1,Y11,Y1p1
// LebedevStencil(5): Y00, Y1m1,Y11,Y1p1, Y2m2,Y2m1,Y20,Y2p1,Y2p2
// LebedevStencil(7): Y00, Y1m1,Y11,Y1p1, Y2m2,Y2m1,Y20,Y2p1,Y2p2, Y3m3,Y3m2,Y3m1,Y30,Y3p1,Y3p2,Y3p3
// nCoefficients = 1+3+5+...+nOrder is the number of spherical harmonics which can be integrated exactly,
// when flattenign the sphercal harmonics: i = l * (l + 1) + m
struct LebedevStencil : public Stencil
{
    LebedevStencil(size_t nOrder);
};



// This type of stencil is not used in the 'general relativistic Lattice Boltzmann Method for radiative transport' code.
// It is only here for reference, if someone might want to test something.
// The 'Stencil' class above is an approximation of this stencil class.
// The approximation allows for a uniform grid which makes interpolations in velocity space easier.
struct GaussLegendreStencil : public Stencil
{
    GaussLegendreStencil(size_t nOrder);
};
#endif //__INCLUDE_GUARD_Stencil_h__