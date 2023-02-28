#ifndef __INCLUDE_GUARD_Radiation_h__
#define __INCLUDE_GUARD_Radiation_h__
#include <math.h>					// Basic math.
#include <omp.h>					// Multithreading.
#include <fstream>					// File input/output.
#include <vector>					// Basic vectors.
#include <string>					// Basic strings.
#include "ControlFlow.hh"			// Template arguments and profiling macros.
#include "Utility.hh"				// Utility functions.
#include "TensorTypes.hh" 			// General relativity tensors.
#include "Grid.h"					// Numerical Grid and mapping to physical domain.
#include "Interpolation.hh"			// Basic interpolation schemes.
#include "Metric.h"				    // Metric parent class.
#include "GeodesicEquationSolver.h"	// Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "TensorOperations.h"       // More specific tensor operations, nullNormalize etc.
#include "Stencil.hh"				// Velocity stencils.
#include "SphericalHarmonics.h"     // Real spherical harmonic functions and expansion.
#include "Log.hh"					// log final results.
#include "Camera.h"					// orthographic camera to take images of radiation field.



// Input system:
enum StreamingType { FlatStatic, FlatDynamic, GeodesicStatic, GeodesicDynamic, CurvedStatic, CurvedDynamic };
std::string StreamingName(int n);



struct Config
{
    std::string name;
    double simTime;
    int writeFrequency;
    bool updateSphericalHarmonics;
    bool keepSourceNodesActive;
    bool writeData;
    bool printToTerminal;
	bool useCamera;
};



class Radiation
{
public:
	StreamingType streamingType;
	Grid& grid;
	Metric& metric;
	Stencil& stencil;
	LebedevStencil& lebedevStencil;
	Camera& camera;


	bool* isInitialGridPoint;
	RealBuffer initialE;
	RealBuffer initialNx;
	RealBuffer initialNy;
	RealBuffer initialNz;
	RealBuffer initialKappa0;
	RealBuffer initialKappa1;
	RealBuffer initialKappaA;
	RealBuffer initialEta;

	RealBuffer nx;
	RealBuffer ny;
	RealBuffer nz;
	RealBuffer nxNew;
	RealBuffer nyNew;
	RealBuffer nzNew;
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
	RealBuffer InorthNew;
	RealBuffer IsouthNew;
	RealBuffer coefficientsS;
	RealBuffer coefficientsX;
	RealBuffer coefficientsY;
	RealBuffer coefficientsZ;
	RealBuffer coefficientsCx;
	RealBuffer coefficientsCy;
	RealBuffer coefficientsCz;

	double normThreshhold = 1e-4;

	Radiation() = delete;
	Radiation(Metric& metric, Stencil& stencil, LebedevStencil& lebedevStencil, Camera& camera, StreamingType streamingType);
	~Radiation();

	int Index(int ijk, int d);
	int Index(int ijk, int d0, int d1);
	int Index(int i, int j, int k, int d);
	int Index(int i, int j, int k, int d0, int d1);

	int HarmonicIndex(int f, int ijk);

	Coord GetTempCoordinate(int ijk, double theta, double phi);
	Tensor3 GetTemp3Velocity(int ijk, double theta, double phi);
	double GetFrequencyShift(int ijk, double theta, double phi);
	double IntensityAt(int ijk, Tensor3 vTempIF);
	Tensor3 AverageF(int i, int j, int k);

	void NormalizeInitialDirections();
	void LoadInitialData();
	void NormalizeInitialIntensities();
	void UpdateSphericalHarmonicsCoefficients();
	void ComputeMomentsIF();
	void ComputeMomentsLF();
	void SetIntensitiesNorthSouth();
	void StreamCurvedStatic();
	void StreamCurvedDynamic();
	void StreamFlatStatic();
	void StreamFlatDynamic();
	void StreamGeodesicStatic();
	void StreamGeodesicDynamic();
	void Collide();
	void TakePicture();

	void RunSimulation(Config config);
};
#endif //__INCLUDE_GUARD_Radiation_h__