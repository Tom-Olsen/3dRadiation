#ifndef __INCLUDE_GUARD_Includes_hh__
#define __INCLUDE_GUARD_Includes_hh__

// Files ending with .h have a corresponding .cpp file and can thus be compiled seperately.
// Files ending with .hh are include only.
// The compilation steps are completely handeld by the CMake/Makefile.
// It is sufficient to include this 'Includes.hh' file in your cpp file to use the code.

#include "Profiler.hh"              // Time measurement profiler.
#include "ControlFlow.hh"           // Template arguments and profiling macros.
#include "Utility.hh"               // Utility functions.
#include "TensorTypes.hh"           // General relativity tensors.
#include "Grid.h"                   // Numerical Grid and mapping to physical domain.
#include "Interpolation.hh"         // Basic interpolation schemes.
#include "Metric.h"                 // Metric parent class.
#include "Spacetimes.h"             // Metric data.
#include "GeodesicEquationSolver.h" // Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "TensorOperations.h"       // More specific tensor operations, nullNormalize etc.
#include "Stencil.h"                // Velocity stencil.
#include "SphericalHarmonics.h"     // Real spherical harmonic functions and expansion.
#include "Radiation.h"              // 'General Relativistic Lattice Boltzmann Method for Radiative Transport' code.
#include "Camera.h"                 // Orthographic Camera to grab 2d slices of intensities orthogonal to camera plane.

#endif //__INCLUDE_GUARD_Includes_hh__