#ifndef __INCLUDE_GUARD_Radiation_h__
#define __INCLUDE_GUARD_Radiation_h__
#include "GeodesicEquationSolver.h" // Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "AdvancedMath.h"       // More specific tensor operations, nullNormalize etc.
#include "SphereGrid.h"             // Helps interpolation on sphere.
#include "SphericalHarmonics.h"     // Real spherical harmonic functions and expansion.
#include "Log.hh"                   // log final results.
#include "Camera.h"                 // orthographic camera to take images of radiation field.
#include "Config.hh"                // Config for simulation parameters.

/*
class Radiation
{
public:
	StreamingType streamingType;
	Grid& grid;
	Metric& metric;
	MyStencil& stencil;
	LebedevStencil& lebedevStencil;
	Camera& camera;
	double sigma = 1.0;


	bool* isInitialGridPoint;
	RealBuffer initialE;
	RealBuffer initialNx;
	RealBuffer initialNy;
	RealBuffer initialNz;
	RealBuffer initialKappa0;
	RealBuffer initialKappa1;
	RealBuffer initialKappaA;
	RealBuffer initialEta;

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
    RealBuffer Inorth;
    RealBuffer Isouth;
	RealBuffer coefficientsS;
	RealBuffer coefficientsX;
	RealBuffer coefficientsY;
	RealBuffer coefficientsZ;
	RealBuffer coefficientsCx;
	RealBuffer coefficientsCy;
	RealBuffer coefficientsCz;

	Radiation() = delete;
	Radiation(Metric& metric, MyStencil& stencil, LebedevStencil& lebedevStencil, Camera& camera, StreamingType streamingType);
	~Radiation();

	size_t Index(size_t ijk, size_t d);
	size_t Index(size_t ijk, size_t d0, size_t d1);
	size_t Index(size_t i, size_t j, size_t k, size_t d);
	size_t Index(size_t i, size_t j, size_t k, size_t d0, size_t d1);

	size_t HarmonicIndex(size_t f, size_t ijk);

	Coord GetTempCoordinate(size_t ijk, Tensor3 direction);
	Tensor3 GetTemp3VelocityIF(size_t ijk, Tensor3 direction);
	double GetFrequencyShift(size_t ijk, Tensor3 direction);
	double IntensityAt(size_t ijk, Tensor3 vTempIF);
	Tensor3 AverageF(size_t i, size_t j, size_t k);

	void NormalizeInitialDirections();
	void LoadInitialData();
	void NormalizeInitialIntensities();
	void UpdateSphericalHarmonicsCoefficients();
	void ComputeMomentsIF();
	void ComputeMomentsLF();
    void SetPoleIntensities();
	void UpdateQuaternions();

	template<class StaticOrDynamic>
	void StreamFlatKernal(size_t i, size_t j, size_t k, size_t d0, size_t d1);
	template<class StaticOrDynamic>
	void StreamCurvedKernal(size_t i, size_t j, size_t k, size_t d0, size_t d1);

	void StreamFlatStatic();
	void StreamFlatDynamic();
	void StreamCurvedStatic();
	void StreamCurvedDynamic();
	void Collide();
	void TakePicture();
	void WriteIntensitiesToCsv(double time, const int frameNumber, std::string directory, std::string name);

	void RunSimulation(Config config);
};
*/
#endif //__INCLUDE_GUARD_Radiation_h__