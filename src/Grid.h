#ifndef __INCLUDE_GUARD_Grid_h__
#define __INCLUDE_GUARD_Grid_h__
#include <iomanip>              // std::setprecision.
#include <fstream>              // File input/output.
#include <jsoncpp/json/json.h>  // Everything about json files.
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
    int nx;
    int ny;
    int nz;
    int nxy;
    int nxyz;
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
    Grid(int nx, int ny, int nz, Coord start, Coord end);
    Grid(const Grid& grid);

    // Setters/Getters:
    void SetCFL(double cfl);
    double GetCFL();

    // Grid Access Tools:
    int Index(int i, int j, int k);
    int Index(int i, int j, int k, int d);
    Coord xyz(int i, int j, int k);
    Coord xyz(double i, double j, double k);
    double i(double x);
    double j(double y);
    double k(double z);
    Coord ijk(const Coord& xyz);

    // Domain Checks:
    bool OutsideDomain(const Coord& xyz);
    bool OutsideDomain(double i, double j, double k);

    // Write Data to file:
    void WriteFrametoJson
    (float time, double* r, double* g, double* b, double* a,
     const int frameNumber, std::string directory, std::string name="");
    void WriteFrametoCsv
    (float time, double* r, double* g, double* b, double* a,
     const int frameNumber, std::string directory, std::string name="");
};
#endif //__INCLUDE_GUARD_Grid_h__