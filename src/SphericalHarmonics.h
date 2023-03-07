#ifndef __INCLUDE_GUARD_SphericalHarmonics_h__
#define __INCLUDE_GUARD_SphericalHarmonics_h__
#include <math.h>           // Basic math.
#include <vector>           // Basic vector.
#include <iostream>         // Output to terminal.
#include "Utility.hh"       // Utility functions.
#include "Stencil.hh"       // Velocity stencil.
#include "TensorTypes.hh"   // General relativity tensors.

// This namespace holds the Spherical Harmonics Ylm
// Note that inputs (x,y,z) or xyz must always be normalized!
class SphericalHarmonics
{
public:
    // l = 0:
    static double Y00(double x, double y, double z);
    static double Y00(const Tensor3& xyz);
    
    // l = 1:
    static double Y1m1(double x, double y, double z);
    static double Y10 (double x, double y, double z);
    static double Y1p1(double x, double y, double z);
    static double Y1m1(const Tensor3& xyz);
    static double Y10 (const Tensor3& xyz);
    static double Y1p1(const Tensor3& xyz);
    
    // l = 2:
    static double Y2m2(double x, double y, double z);
    static double Y2m1(double x, double y, double z);
    static double Y20 (double x, double y, double z);
    static double Y2p1(double x, double y, double z);
    static double Y2p2(double x, double y, double z);
    static double Y2m2(const Tensor3& xyz);
    static double Y2m1(const Tensor3& xyz);
    static double Y20 (const Tensor3& xyz);
    static double Y2p1(const Tensor3& xyz);
    static double Y2p2(const Tensor3& xyz);
    
    // l = 3:
    static double Y3m3(double x, double y, double z);
    static double Y3m2(double x, double y, double z);
    static double Y3m1(double x, double y, double z);
    static double Y30 (double x, double y, double z);
    static double Y3p1(double x, double y, double z);
    static double Y3p2(double x, double y, double z);
    static double Y3p3(double x, double y, double z);
    static double Y3m3(const Tensor3& xyz);
    static double Y3m2(const Tensor3& xyz);
    static double Y3m1(const Tensor3& xyz);
    static double Y30 (const Tensor3& xyz);
    static double Y3p1(const Tensor3& xyz);
    static double Y3p2(const Tensor3& xyz);
    static double Y3p3(const Tensor3& xyz);
    
    // l = 4:
    static double Y4m4(double x, double y, double z);
    static double Y4m3(double x, double y, double z);
    static double Y4m2(double x, double y, double z);
    static double Y4m1(double x, double y, double z);
    static double Y40 (double x, double y, double z);
    static double Y4p1(double x, double y, double z);
    static double Y4p2(double x, double y, double z);
    static double Y4p3(double x, double y, double z);
    static double Y4p4(double x, double y, double z);
    static double Y4m4(const Tensor3& xyz);
    static double Y4m3(const Tensor3& xyz);
    static double Y4m2(const Tensor3& xyz);
    static double Y4m1(const Tensor3& xyz);
    static double Y40 (const Tensor3& xyz);
    static double Y4p1(const Tensor3& xyz);
    static double Y4p2(const Tensor3& xyz);
    static double Y4p3(const Tensor3& xyz);
    static double Y4p4(const Tensor3& xyz);

    static std::vector<double (*)(double, double, double)> Harmonics;

    static std::vector<double> GetCoefficients(const Stencil& stencil, const double* data, int nCoefficients);
    static std::vector<double> GetCoefficients(const LebedevStencil& stencil, const double* data, int nCoefficients);
    static std::vector<double> GetCoefficients(const GaussLegendreStencil& stencil, const double* data, int nCoefficients);
    static void GetCoefficients(const Stencil& stencil, const double* data, int nCoefficients, double* coefficients);
    static void GetCoefficients(const LebedevStencil& stencil, const double* data, int nCoefficients, double* coefficients);
    static void GetCoefficients(const GaussLegendreStencil& stencil, const double* data, int nCoefficients, double* coefficients);
    static double GetValue(Tensor3 direction, const std::vector<double>& coefficients, int nCoefficients);
    static double GetValue(Tensor3 direction, double* coefficients, int nCoefficients);
    static double GetValue(double theta, double phi, const std::vector<double>& coefficients, int nCoefficients);
    static double GetValue(double theta, double phi, double* coefficients, int nCoefficients);
};
#endif //__INCLUDE_GUARD_SphericalHarmonics_h__