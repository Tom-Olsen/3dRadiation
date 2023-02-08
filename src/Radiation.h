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
#include "Log.hh"					// log final results



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
};



class Radiation
{
public:
	StreamingType streamingType;
	Grid& grid;
	Metric& metric;
	Stencil& stencil;
	LebedevStencil& lebedevStencil;

	bool* isInitialGridPoint;
	double* initialE;
	double* initialNx;
	double* initialNy;
	double* initialNz;
	double* initialKappa0;
	double* initialKappa1;
	double* initialKappaA;
	double* initialEta;

	double* nx;
	double* ny;
	double* nz;
	double* nxNew;
	double* nyNew;
	double* nzNew;
	double* E;
	double* Fx;
	double* Fy;
	double* Fz;
	double* Pxx;
	double* Pxy;
	double* Pxz;
	double* Pyy;
	double* Pyz;
	double* Pzz;
	double* E_LF;
	double* Fx_LF;
	double* Fy_LF;
	double* Fz_LF;
	double* Pxx_LF;
	double* Pxy_LF;
	double* Pxz_LF;
	double* Pyy_LF;
	double* Pyz_LF;
	double* Pzz_LF;
    double* kappa0;
    double* kappa1;
    double* kappaA;
    double* eta;
	double* I;
	double* Inew;
	double* Inorth;
	double* Isouth;
	double* coefficientsS;
	double* coefficientsX;
	double* coefficientsY;
	double* coefficientsZ;
	double* coefficientsCx;
	double* coefficientsCy;
	double* coefficientsCz;

	double normThreshhold = 1e-4;

	Radiation() = delete;
	Radiation(Metric& metric_, Stencil& stencil_, LebedevStencil& lebedevStencil_, StreamingType streamingType_);
	~Radiation();

	Coord GetTempCoordinate(int i, int j, int k, double theta, double phi);
	Tensor3 GetTemp3Velocity(int i, int j, int k, double theta, double phi);
	double GetFrequencyShift(int i, int j, int k, double theta, double phi);
	double IntensityAt(int ijk, Tensor3 vTempIF);
	Tensor3 AverageF(int i, int j, int k);

	void NormalizeInitialDirections();
	void LoadInitialData();
	void NormalizeInitialData();
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

	void RunSimulation(Config config);
};
#endif //__INCLUDE_GUARD_Radiation_h__