#ifndef __INCLUDE_GUARD_Radiation_h__
#define __INCLUDE_GUARD_Radiation_h__
#include "GeodesicEquationSolver.h" // Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "SpecialMath.h"            // More specific tensor operations, nullNormalize etc.
#include "SphericalHarmonics.h"     // Real spherical harmonic functions and expansion.
#include "InterpolationGrid.h"      // Speed up unstructured grid interpolation in velocity space.
#include "Logger.hh"                // logs final results.
#include "Camera.h"                 // orthographic camera to take images of radiation field.
#include "Config.hh"                // Config for simulation parameters.

class Radiation
{
private:
    // Constants:
    static constexpr int HALO = 1;
    static constexpr double MIN_FLUX_NORM = 1e-16;
    static constexpr double MIN_ENERGY_DENSITY = 1e-16;
    static constexpr double LAMBDA_ITTERATION_TOLERENCE = 1e-12;
    static constexpr int MAX_LAMBDA_ITERATIONS = 100;
    static constexpr double MAX_INTERPOLATION_ERROR = 0.20; // in %
    static constexpr glm::vec3 from = glm::vec3(0, 0, 1);
    // Units conversion is not correct yet!
    double etaCGStoCode = 7.67822;
    double kappaCGStoCode = 1.47760;

public:
    Grid &grid;
    Metric &metric;
    Stencil &stencil;
    LebedevStencil &streamingStencil;
    InterpolationGrid &interpGrid;
    Camera &camera;
    Config config;
    Logger logger;

    // Set from the outside:
    bool *isInitialGridPoint;
    RealBuffer initialE_LF;
    RealBuffer initialFx_LF;
    RealBuffer initialFy_LF;
    RealBuffer initialFz_LF;
    RealBuffer initialPxx_LF;
    RealBuffer initialPxy_LF;
    RealBuffer initialPxz_LF;
    RealBuffer initialPyy_LF;
    RealBuffer initialPyz_LF;
    RealBuffer initialPzz_LF;
    RealBuffer initialKappa0;
    RealBuffer initialKappa1;
    RealBuffer initialKappaA;
    RealBuffer initialEta;
    RealBuffer initialI;
    QuatBuffer initialQ;
    // kappa* and eta must be givne in CGS units
    // and will be converted to code units in the LoadInitialData() method.

    RealBuffer sigma;
    RealBuffer normalization;
    QuatBuffer q;
    QuatBuffer qNew;
    RealBuffer E;
    RealBuffer Fx;
    RealBuffer Fy;
    RealBuffer Fz;
    RealBuffer Pxx;
    RealBuffer Pxy;
    RealBuffer Pxz;
    RealBuffer Pyy;
    RealBuffer Pyz;
    RealBuffer Pzz;
    RealBuffer E_LF;
    RealBuffer Fx_LF;
    RealBuffer Fy_LF;
    RealBuffer Fz_LF;
    RealBuffer Pxx_LF;
    RealBuffer Pxy_LF;
    RealBuffer Pxz_LF;
    RealBuffer Pyy_LF;
    RealBuffer Pyz_LF;
    RealBuffer Pzz_LF;
    RealBuffer kappa0;
    RealBuffer kappa1;
    RealBuffer kappaA;
    RealBuffer eta;
    RealBuffer I;
    RealBuffer Inew;
    RealBuffer coefficientsS;
    RealBuffer coefficientsX;
    RealBuffer coefficientsY;
    RealBuffer coefficientsZ;
    RealBuffer coefficientsCx;
    RealBuffer coefficientsCy;
    RealBuffer coefficientsCz;

    Radiation() = delete;
    Radiation(Metric &metric, Stencil &stencil, LebedevStencil &streamingStencil, InterpolationGrid &interpGrid, Camera &camera, Config config);
    ~Radiation();

    size_t Index(size_t ijk, size_t d);
    size_t Index(size_t i, size_t j, size_t k, size_t d);

    size_t HarmonicIndex(size_t f, size_t ijk);

    Coord GetTempCoordinate(size_t ijk, const Tensor3 &direction);
    Tensor3 GetTemp3VelocityIF(size_t ijk, const Tensor3 &direction);
    double GetFrequencyShift(size_t ijk, const Tensor3 &direction);
    double IntensityAt(size_t ijk, Tensor3 vTempIF);
    Tensor3 AverageF(size_t i, size_t j, size_t k);

    Tensor4 InitialDataLFtoIF(size_t ijk);
    void InitSigmaAndNormalization();
    double SigmaMax();
    double FluxMax();
    void LoadInitialData();
    void UpdateSphericalHarmonicsCoefficients();
    void ComputeMomentsIF();
    void ComputeMomentsLF();
    void UpdateQuaternions();

    void StreamFlatFixed();
    void StreamFlatAdaptive();
    void StreamCurvedFixed();
    void StreamCurvedAdaptive();

    void CollideStaticFluidForwardEuler();
    void CollideStaticFluidBackwardEuler();
    void CollideForwardEuler();
    void CollideBackwardEuler();

    void TakePicture();
    void RunSimulation();
};
#endif //__INCLUDE_GUARD_Radiation_h__