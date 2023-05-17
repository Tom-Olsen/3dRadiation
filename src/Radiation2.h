#ifndef __INCLUDE_GUARD_Radiation2_h__
#define __INCLUDE_GUARD_Radiation2_h__
#include "GeodesicEquationSolver.h" // Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "AdvancedMath.h"       // More specific tensor operations, nullNormalize etc.
#include "SphereGrid.h"             // Helps interpolation on sphere.
#include "SphericalHarmonics.h"     // Real spherical harmonic functions and expansion.
#include "Log.hh"                   // log final results.
#include "Camera.h"                 // orthographic camera to take images of radiation field.
#include "Config.hh"                // Config for simulation parameters.

// This version uses spherical harmonic coefficiens instead of intensity distributions.
// I[ijk,d] = SphericalHarmonicsXyz::GetValue(C(d), &coefficientsI[ijk], nCoefficients)
// Is the actual intensity I(x,c) at position x, in direction c.
// It is not an intesity population I_d(x) = w(d) * I(x,c)
// => E = sum_d w_d I(x,c_d) instead of E = sum_d I_d
// Sadly it is way to slow.

class Radiation2
{
public:
	Grid& grid;
	Metric& metric;
	Stencil& intensityStencil;
	Stencil& streamingStencil;
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
	RealBuffer coefficientsI;
	RealBuffer coefficientsInew;
	RealBuffer coefficientsS;
	RealBuffer coefficientsX;
	RealBuffer coefficientsY;
	RealBuffer coefficientsZ;
	RealBuffer coefficientsCx;
	RealBuffer coefficientsCy;
	RealBuffer coefficientsCz;

	Radiation2() = delete;
	Radiation2(Metric& metric, Stencil& intensityStencil, Stencil& streamingStencil, Camera& camera);
	~Radiation2();

	size_t IntensityHarmonicIndex(size_t ijk, size_t d);
	size_t IntensityHarmonicIndex(size_t i, size_t j, size_t k, size_t d);
	size_t StreamingHarmonicIndex(size_t ijk, size_t d);

	void NormalizeInitialDirections();
	void LoadInitialData();
	void NormalizeInitialIntensities();
	void UpdateStreamingCoefficients();
	void ComputeMomentsIF();
	void ComputeMomentsLF();

	Coord GetTempCoordinate(size_t ijk, Tensor3 direction);
	Tensor3 GetTemp3VelocityIF(size_t ijk, Tensor3 direction);
	double GetFrequencyShift(size_t ijk, Tensor3 direction);
	double IntensityAt(size_t ijk, const Tensor3& vTempIF);

	void Stream();
	void Collide();
	void TakePicture();

	void RunSimulation(Config config);
};
#endif //__INCLUDE_GUARD_Radiation2_h__