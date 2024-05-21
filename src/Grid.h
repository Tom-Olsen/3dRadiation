#ifndef __INCLUDE_GUARD_Grid_h__
#define __INCLUDE_GUARD_Grid_h__
#include <fstream>        // File input/output.
#include "ControlFlow.hh" // Template arguments and profiling macros.
#include "Utility.hh"     // Utility functions.
#include "DataTypes.hh"   // General relativity tensors.
#include "Profiler.hh"    // Time measurement profiler.

class Grid
{
private:
    double m_cfl = 1;

public:
    // Members:
    size_t nx, ny, nz, nxy, nxyz;
    double dx, dy, dz, dt;
    double startx, starty, startz;
    double endx, endy, endz;

    // Constructors:
    Grid() = delete;
    Grid(size_t nx, size_t ny, size_t nz, Coord start, Coord end);
    Grid(const Grid &grid);

    // Setters/Getters:
    void SetCFL(double cfl);
    double GetCFL();

    // Grid Access Tools:
    size_t Index(size_t i, size_t j, size_t k);
    double x(size_t i);
    double y(size_t j);
    double z(size_t k);
    double x(double i);
    double y(double j);
    double z(double k);
    Coord xyz(size_t i, size_t j, size_t k);
    Coord xyz(double i, double j, double k);
    Coord xyz(size_t ijk);
    double i(double x);
    double j(double y);
    double k(double z);
    size_t i(size_t ijk);
    size_t j(size_t ijk);
    size_t k(size_t ijk);
    Coord ijk(const Coord &xyz);

    // Domain Checks:
    bool OutsideDomain(const Coord &xyz);
    bool OutsideDomain(double i, double j, double k);

    // Write Data to file:
    void WriteFrametoCsv(float time, const RealBuffer &r, const RealBuffer &g, const RealBuffer &b, const RealBuffer &a, std::string directory, std::string name = "");

    // Debugging:
    void Print();
};
#endif //__INCLUDE_GUARD_Grid_h__