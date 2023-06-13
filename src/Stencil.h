#ifndef __INCLUDE_GUARD_Stencil_h__
#define __INCLUDE_GUARD_Stencil_h__
#include <vector>
#include <set>
#include <span>
#include "Utility.hh"
#include "DataTypes.hh"
#include "ConvexHull.h"
#include "Interpolation.h"
#include "SpecialMath.h"



struct Stencil;



struct InterpolationGrid
{
    // Public Properties:
    size_t nTh;
    size_t nPh;
    size_t nDir;

    // Constructor:
    InterpolationGrid() = default;
    InterpolationGrid(size_t nTh, size_t nPh);
    
    // Grid operations:
    double Theta(size_t i) const;
    double Phi(size_t j) const;
    double i(double theta) const;
    double j(double phi) const;
    size_t Index(size_t i, size_t j) const;
    Tensor3 C(size_t i, size_t j) const;
};



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
    std::string name;       // Name of the stencel, e.g. "Lebedev9.1.4"
    size_t nDir;            // Number of directions in stencil.
    size_t nOrder;          // Quadrature integration order.
    size_t nCoefficients;   // Number of exactly integrated Spherical Harmonics, counted with flat index i = l * (l + 1) + m.
    size_t nRings;          // Number of ghost rings to be added.
    size_t nRing0;          // Number of ghost directions in smallest ring.
    size_t nGhost;          // Total amount of ghost directions.
    double thetaGhost;      // Ghost rings are evenly spread from in (0,thetaGhost), excluding boundaries.
    Mesh mesh;
    InterpolationGrid interpolationGrid;
private:
    int interpolationGridRes = 1000;

protected: // Internal Buffers:
    RealBuffer w;
    RealBuffer cx;
    RealBuffer cy;
    RealBuffer cz;
    RealBuffer theta;
    RealBuffer phi;
    
    UnstructuredMatrix<Vector3Int> connectedTriangles;
    UnstructuredMatrix<Vector3> voronoiCells;
    UnstructuredMatrix<size_t> voronoiNeighbours;

protected: // Internal Buffers on interpolationGrid:
    SizeTBuffer neighbour0OnGrid;
    SizeTBuffer neighbour1OnGrid;
    SizeTBuffer neighbour2OnGrid;
public: // Internal Buffers on interpolationGrid:
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
    void AllocateBuffers();
    void AddGhostDirections();
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
    LebedevStencil(size_t nOrder, size_t nRings=0, size_t nRing0=0, double thetaGhost=0);
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