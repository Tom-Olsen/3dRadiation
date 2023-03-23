#ifndef __INCLUDE_GUARD_SphericalHarmonics_h__
#define __INCLUDE_GUARD_SphericalHarmonics_h__
#include <math.h>           // Basic math.
#include <vector>           // Basic vector.
#include <iostream>         // Output to terminal.
#include "Utility.hh"       // Utility functions.
#include "Stencil.h"        // Velocity stencil.
#include "TensorTypes.hh"   // General relativity tensors.



struct SphericalHarmonicsThPh
{
private:
    // Associated Legendre Polynomials.
    template<int L, int M>
    static double P(double x);

    // Renormalisation constant for Spherical Harmonic function.
    template<int L, int M>
    static double Kd();

    // Harmonics Templates:
    template<int L, int M>
    static double Y(double theta, double phi);
    template<int I>
    static double Y(double theta, double phi);
    
    // Function Pointer:
    static std::vector<double (*)(double, double)> Ylist;
public:
    static double Y(int i, double theta, double phi);
    static double Y(int l, int m, double theta, double phi);
    

    // Spherical Harmonics Expansion:
    static std::vector<double> GetCoefficients(const Stencil& stencil, const double* data);
    static void GetCoefficients(const Stencil& stencil, const double* data, double* coefficients);
    static double GetValue(double theta, double phi, const std::vector<double>& coefficients, size_t nCoefficients);
    static double GetValue(double theta, double phi, double* coefficients, size_t nCoefficients);
};



struct SphericalHarmonicsXyz
{
private:
    // Harmonics Template:
    template<int I>
    static double Y(double x, double y, double z);

    // Function Pointer:
    static std::vector<double (*)(double, double, double)> Ylist;
public:
    static double Y(int i, double x, double y, double z);
    static double Y(int l, int m, double x, double y, double z);
    static double Y(int i, const Tensor3& dir);
    static double Y(int l, int m, const Tensor3& dir);
    
    // Spherical Harmonics Expansion:
    static std::vector<double> GetCoefficients(const Stencil& stencil, const double* data);
    static void GetCoefficients(const Stencil& stencil, const double* data, double* coefficients);
    static double GetValue(double x, double y, double z, const std::vector<double>& coefficients, size_t nCoefficients);
    static double GetValue(double x, double y, double z, double* coefficients, size_t nCoefficients);
    static double GetValue(const Tensor3& direction, const std::vector<double>& coefficients, size_t nCoefficients);
    static double GetValue(const Tensor3& direction, double* coefficients, size_t nCoefficients);
};
#endif //__INCLUDE_GUARD_SphericalHarmonics_h__