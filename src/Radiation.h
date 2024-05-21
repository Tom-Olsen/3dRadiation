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
    static constexpr double LAMBDA_ITTERATION_TOLERENCE = 1e-4; // 0.01% error
    static constexpr int MAX_LAMBDA_ITERATIONS = 100;
    static constexpr glm::vec3 from = glm::vec3(0, 0, 1);

public:
    Grid &grid;
    Metric &metric;
    LebedevStencil &stencil;
    LebedevStencil &streamingStencil;
    InterpolationGrid &interpGrid;
    Camera &camera;
    Config config;
    Logger logger;

    // Initial data is set from the outside (in code units):
    bool *isInitialGridPoint;
    // config.initialDataType = InitialDataType::Intensities
    RealBuffer initialI;
    // config.initialDataType = InitialDataType::Moments
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

    // Rotation of stencils:
    QuatBuffer q;
    QuatBuffer qNew;

    // Inertial frame moments:
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

    // Lab frame moments:
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

    // Fluid properties:
    RealBuffer kappa0;
    RealBuffer kappa1;
    RealBuffer kappaA;
    RealBuffer eta;
    RealBuffer ux;
    RealBuffer uy;
    RealBuffer uz;

    // Population intensities:
    RealBuffer I;
    RealBuffer Inew;
    
    // Spherical harmonics coefficients for geodesic streaming:
    RealBuffer coefficientsS;
    RealBuffer coefficientsX;
    RealBuffer coefficientsY;
    RealBuffer coefficientsZ;
    RealBuffer coefficientsCx;
    RealBuffer coefficientsCy;
    RealBuffer coefficientsCz;

    // Testing:
    IntBuffer itterationCount;
    int maxItterationCount = 0;
    double averageItterationCount = 0;

    Radiation() = delete;
    Radiation(Metric &metric, LebedevStencil &stencil, LebedevStencil &streamingStencil, InterpolationGrid &interpGrid, Camera &camera, Config config);
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
    void LoadInitialData();
    void LoadInitialDataIntensitiesFixed();
    void LoadInitialDataMomentsFixed();
    void LoadInitialDataIntensitiesAdaptive();
    void LoadInitialDataMomentsAdaptive();
    void UpdateSphericalHarmonicsCoefficients();
    void ComputeMomentsIF();
    void ComputeMomentsLF();
    void UpdateQuaternions();

    void StreamFlatFixed();
    void StreamFlatAdaptive();
    void StreamCurvedFixed();
    void StreamCurvedAdaptive();

    void Collide();
    //void CollideStaticFluidForwardEuler();
    //void CollideStaticFluidBackwardEuler();
    //void CollideForwardEuler();
    //void CollideBackwardEuler();

    void TakePicture();
    void RunSimulation();
};
#endif //__INCLUDE_GUARD_Radiation_h__