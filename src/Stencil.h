#ifndef __INCLUDE_GUARD_Stencil_h__
#define __INCLUDE_GUARD_Stencil_h__
#include <vector>
#include <set>
#include <span>
#include <tuple>
#include <fstream>
#include "Utility.hh"
#include "DataTypes.hh"
#include "ConvexHull.h"
#include "Interpolation.h"
#include "SpecialMath.h"

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
public:
    // Public Properties:
    std::string name;            // Name of the stencel, e.g. "Lebedev9.1.4"
    size_t nDir;                 // Number of directions in stencil.
    size_t nOrder;               // Quadrature integration order.
    size_t nCoefficients;        // Number of exactly integrated Spherical Harmonics, counted with flat index i = l * (l + 1) + m.
    size_t nGhost;               // Number of ghost directions.
    double refinement0Threshold; // Number of ghost rings to be added.
    double refinement1Threshold; // Number of ghost directions in smallest ring.
    double refinement2Threshold; // Total amount of ghost directions.
    LookUpTable sigmaOfRelativeFlux;
    double sigmaMax;
    double relativeFluxMax;
    Mesh mesh;
    static constexpr double maxInterpolationError  = 0.01; // 1% (standard value)
    // static constexpr double maxInterpolationError  = 0.20; // (test higher values)
    static constexpr double sigmaMaxSearchAccuracy = 1e-3;

protected:
    // Internal Buffers (data for each vertex):
    RealBuffer w;
    RealBuffer cx;
    RealBuffer cy;
    RealBuffer cz;
    RealBuffer theta;
    RealBuffer phi;
    UnstructuredMatrix<Vector3Int> connectedTriangles;
    UnstructuredMatrix<Vector3> voronoiCells;
    UnstructuredMatrix<size_t> voronoiNeighbours;

public:
    // Getters:
    double W(size_t d) const;
    double Cx(size_t d) const;
    double Cy(size_t d) const;
    double Cz(size_t d) const;
    double Theta(size_t d) const;
    double Phi(size_t d) const;
    Tensor3 Ct3(size_t d) const;
    Vector3 Cv3(size_t d) const;
    std::span<const Vector3Int> ConnectedTriangles(size_t d) const;
    std::span<const Vector3> VoronoiCell(size_t d) const;
    std::span<const size_t> VoronoiNeighbours(size_t d) const;

    // Expensive Getters:
    // Get voronoi neighbours and weights of an abitrary point p.
    std::tuple<std::vector<size_t>, std::vector<double>> VoronoiNeighboursAndWeights(const Vector3 &p) const;
    std::tuple<std::vector<size_t>, std::vector<double>> VoronoiNeighboursAndWeights(const Vector3 &p, bool test) const;
    // Get barycentric neighbours and weights of an abitrary point p.
    std::tuple<Vector3Int, Vector3> BarycentricNeighboursAndWeights(const Vector3 &p) const;

private:
    // Internal helper methods:
    // Index of nearest stencil.C(?) for arbitrary point p.
    size_t NearestNeighbour(const Vector3 &p) const;
    // Creates virtual voronoi cell of arbitrary point p. This includes the natural neighbours and voronoi points that are inside the virual cell.
    std::vector<Vector3> VirtualVoronoiCellOf(const Vector3 &p, std::vector<size_t> &naturalNeighbours, std::vector<Vector3> &pointsInsideCell) const;
    // Get voronoi weights of virtual voronoi cell.
    std::vector<double> VoronoiWeights(const std::vector<Vector3> &virtualCell, const std::vector<size_t> &naturalNeighbours, const std::vector<Vector3> &pointsInsideCell) const;

protected:
    // Initialization:
    void AllocateBuffers();
    void AddGhostDirections();
    void SortDirections();
    void InitializeMesh();
    void InitializeConnectedTriangles();
    void InitializeVoronoiCells();
    void InitializeVoronoiNeighbours();
    void PopulateLookUpTable();
    double MaxSigma();

public:
    // Debugging:
    void Print() const;
    void PrintAll() const;
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
    LebedevStencil(size_t nOrder, double refinement0Threshold = 0, double refinement1Threshold = 0, double refinement2Threshold = 0);
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