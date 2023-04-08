#ifndef __INCLUDE_GUARD_Grid_h__
#define __INCLUDE_GUARD_Grid_h__
// #include <iomanip>              // std::setprecision.
#include <fstream>              // File input/output.
#include "ControlFlow.hh"       // Template arguments and profiling macros.
#include "Utility.hh"           // Utility functions.
#include "TensorTypes.hh"       // General relativity tensors.
#include "Profiler.hh"          // Time measurement profiler.



class Grid
{
private:
    double m_cfl=1;
public:
    // Members:
    size_t nx;
    size_t ny;
    size_t nz;
    size_t nxy;
    size_t nxyz;
    double dx;
    double dy;
    double dz;
    double dt;
    double startx;
    double starty;
    double startz;
    double endx;
    double endy;
    double endz;

    // Constructors:
    Grid() = delete;
    Grid(size_t nx, size_t ny, size_t nz, Coord start, Coord end);
    Grid(const Grid& grid);

    // Setters/Getters:
    void SetCFL(double cfl);
    double GetCFL();

    // Grid Access Tools:
    size_t Index(size_t i, size_t j, size_t k);
    Coord xyz(size_t i, size_t j, size_t k);
    Coord xyz(double i, double j, double k);
    double i(double x);
    double j(double y);
    double k(double z);
    size_t i(size_t ijk);
    size_t j(size_t ijk);
    size_t k(size_t ijk);
    Coord ijk(const Coord& xyz);

    // Domain Checks:
    bool OutsideDomain(const Coord& xyz);
    bool OutsideDomain(double i, double j, double k);

    // Write Data to file:
    void WriteFrametoCsv
    (float time, const RealBuffer& r, const RealBuffer& g, const RealBuffer& b, const RealBuffer& a,
     const int frameNumber, std::string directory, std::string name="");
};
#endif //__INCLUDE_GUARD_Grid_h__