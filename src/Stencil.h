#ifndef __INCLUDE_GUARD_Stencil_h__
#define __INCLUDE_GUARD_Stencil_h__
#include <vector>
#include <set>
#include "Utility.hh"
#include "TensorTypes.hh"
#include "ConvexHull.h"
#include "SphereGrid.h"



// Spherical Coordinates convention:
// theta € [0,pi]
//   phi € [0,2pi]
// theta =  0 => North Pole
// theta = pi => South Pole
// (th,phi) = (pi/2,    0) => (x,y,z) = ( 1, 0, 0)
// (th,phi) = (pi/2, pi/2) => (x,y,z) = ( 0, 1, 0)
// (th,phi) = (pi/2,   pi) => (x,y,z) = (-1, 1, 0)
// (th,phi) = (pi/2,3pi/2) => (x,y,z) = ( 0,-1, 0)
// (th,phi) = (pi/2,  2pi) => (x,y,z) = ( 1, 0, 0)



// Base Stencil class. This functionality is supported by all stencils.
struct Stencil
{
public:
    std::string name;
    size_t nDir;                        // Number of directions in stencil.
    size_t nOrder;                      // Quadrature integration.
    size_t nCoefficients;               // Number of exact Spherical Harmonics, counted with flat index i = l * (l + 1) + m.
    // List of lists of all triangles that are connected to one direction vector c.
    UnstructuredMatrix<Vector3Int> connectedTriangles;
    UnstructuredMatrix<size_t> connectedVerticesOrder1;
    UnstructuredMatrix<size_t> connectedVerticesOrder2;
    SphereGrid sphereGrid;
protected:
    RealBuffer w;
    RealBuffer cx;
    RealBuffer cy;
    RealBuffer cz;
    RealBuffer theta;
    RealBuffer phi;
    SizeTBuffer neighbour0;
    SizeTBuffer neighbour1;
    SizeTBuffer neighbour2;
protected:
    void SetCoefficientCount();
    void AllocateBuffers();
    void SortDirections();
    void InitializeConnectedTriangles();
    void InitializeConnectedVertices();
    void InitializeNearestNeighbourGrid();
public:
    virtual double W(size_t d) const;
    virtual double Theta(size_t d) const;
    virtual double Phi(size_t d) const;
    virtual double Cx(size_t d) const;
    virtual double Cy(size_t d) const;
    virtual double Cz(size_t d) const;
    virtual Tensor3 C(size_t d) const;

    size_t GetNeighbourIndex0(const Tensor3& p);
    size_t GetNeighbourIndex1(const Tensor3& p);
    size_t GetNeighbourIndex2(const Tensor3& p);
    Tensor3 GetNeighbour0(const Tensor3& p);
    Tensor3 GetNeighbour1(const Tensor3& p);
    Tensor3 GetNeighbour2(const Tensor3& p);

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