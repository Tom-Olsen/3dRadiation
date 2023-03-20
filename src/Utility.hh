#ifndef __INCLUDE_GUARD_Utility_hh__
#define __INCLUDE_GUARD_Utility_hh__
#include <iomanip>                  // std::setprecision(), std::put:time
#include <iostream>                 // Output to terminal.
#include <fstream>                  // File input/output.
#include <cmath>                    // signbit.
#include <algorithm>                // clamp.
#include <sstream>                  // stringstream.
#include <filesystem>               // Folder/file management.
#include <vector>                   // basic c++ vector
#include "glm/glm/gtc/quaternion.hpp"   // Quaternions.
#include "ControlFlow.hh"           // Template arguments and profiling macros.
#include "eigen/Eigen/Dense"        // Eigen library for solving linear systems.



// Exit programm with Error Message.
inline void exit_on_error(std::string msg="")
{
    std::cout << "ERROR: " + msg + "\n";
    exit(errno);
}



// Memory aligned data buffers:
template <class T>
class AlignedArrayAllocator
{
public:
    using value_type = T;

    // Constructor:
    T* allocate(size_t size)
    { return static_cast<T *>(std::aligned_alloc(64,size * sizeof(T))); }

    // Descturctor:
    void deallocate(T * p_t, size_t)
    { std::free(p_t); }

    // Initialize each element:
    template<class U, class... Args>
    void construct(U* p, Args&&... args)
    {
        if constexpr(std::is_same<U,glm::quat>::value)
            *p = glm::quat(1,0,0,0);
        else
            *p = 0;
    }

    // Needed for std::swap(..., ...):
    bool operator==(const AlignedArrayAllocator &)
    { return true; }
};
using RealBuffer = std::vector<double,AlignedArrayAllocator<double>>;
using QuatBuffer = std::vector<glm::quat,AlignedArrayAllocator<glm::quat>>;



// Fast integer exponentiation:
/// @brief a^N with double a and integer N.
/// @tparam N integer power.
/// @param a number to exponentiate.
/// @return a*a*a...*a, N times.
template<int N>
inline double IntegerPow(double a)
{ return a * IntegerPow<N-1>(a); }
template<>
inline double IntegerPow<0>(double a)
{ return 1; }



/// @brief Polynomial approximation of atan(z). Much faster then std libraries.
/// @tparam Order Polynomial approximation order [0,3,5,7,9,11]. Use 0 for rational approximation.
/// @param z Bound to [-1,1].
/// @return Polynomial approximation of: atan(z)
template<int Order>
inline double MyAtan(double z)
{
    // Rational approximation:
    // maxError = 0.0048829, averageError = 0.0046227
    if constexpr(Order == 0)
    {
        const double a = 1.0;
        const double b = 0.28;
        return z / (a + b * z * z);
    }

    // Polynomial approximations can be found in:
    // Approximations for digital computers,
    // https://blasingame.engr.tamu.edu/z_zCourse_Archive/P620_18C/P620_zReference/PDF_Txt_Hst_Apr_Cmp_(1955).pdf

    // maxError = 0.0049493, averageError = 0.0062662
    if constexpr(Order == 3)
    {
        const double a1 =  0.97239411;
        const double a3 = -0.19194795;
        double zz = z * z;
        return z * (a1 + zz * a3);
    }
    
    // maxError = ???, averageError = ???
    if constexpr(Order == 5)
    {
        const double a1 =  0.995354;
        const double a3 = -0.288679;
        const double a5 =  0.079331;
        double zz = z * z;
        return z * (a1 + zz * (a3 + zz * a5));
    }
    
    // maxError = ???, averageError = ???
    if constexpr(Order == 7)
    {
        const double a1 =  0.9992150;
        const double a3 = -0.3211819;
        const double a5 =  0.1462766;
        const double a7 = -0.0389929;
        double zz = z * z;
        return z * (a1 + zz * (a3 + zz * (a5 + zz * a7)));
    }
    
    // maxError = ???, averageError = ???
    if constexpr(Order == 9)
    {
        const double a1 =  0.9998660;
        const double a3 = -0.3302995;
        const double a5 =  0.1801410;
        const double a7 = -0.0851330;
        const double a9 =  0.0208351;
        double zz = z * z;
        return z * (a1 + zz * (a3 + zz * (a5 + zz * (a7 + zz * a9))));
    }


    // maxError = ???, averageError = ???
    if constexpr(Order == 11)
    {
        const double a1  =  0.99997726;
        const double a3  = -0.33262347;
        const double a5  =  0.19354346;
        const double a7  = -0.11643287;
        const double a9  =  0.05265332;
        const double a11 = -0.01172120;
        double zz = z * z;
        return z * (a1 + zz * (a3 + zz * (a5 + zz * (a7 + zz * (a9 + zz * a11)))));
    }
}
/// @brief Polynomial order given by the macro APPROXIMATION_ORDER in ControlFlow.hh
/// @param z Bound to [-1,1].
/// @return Polynomial approximation of: atan(z)
inline double MyAtan(double x)
{ return MyAtan<APPROXIMATION_ORDER>(x); }



/// @brief Polynomial approximation of atan2(y,x). Much faster then std libraries.
/// @tparam Order Polynomial approximation order [0,3,5,7,9,11].
/// @param y Coordinate on 2D Cartesian plane.
/// @param x Coordinate on 2D Cartesian plane.
/// @return Polynomial approximation of: atan2(y,x)
template<int Order>
inline double MyAtan2(double y, double x)
{// https://www.dsprelated.com/showarticle/1052.php
    constexpr double pi = M_PI;
    constexpr double pi_2 = M_PI_2;
    if (x != 0.0)
    {
        if (fabsf(x) > fabsf(y))
        {
            const double z = y / x;
            if (x > 0.0)
                return MyAtan<Order>(z);
            else if (y >= 0.0)
                return MyAtan<Order>(z) + pi;
            else
                return MyAtan<Order>(z) - pi;
        }
        else // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
        {
            const double z = x / y;
            if (y > 0.0)
                return -MyAtan<Order>(z) + pi_2;
            else
                return -MyAtan<Order>(z) - pi_2;
        }
    }
    else
    {
        if (y > 0.0)
            return pi_2;
        else if (y < 0.0)
            return -pi_2;
    }
    return 0.0; // x,y = 0. Could return NaN instead.
}
/// @brief Polynomial order given by the macro APPROXIMATION_ORDER in ControlFlow.hh
/// @param y Coordinate on 2D Cartesian plane.
/// @param x Coordinate on 2D Cartesian plane.
/// @return Polynomial approximation of: atan2(y,x)
inline double MyAtan2(double y, double x)
{ return MyAtan2<APPROXIMATION_ORDER>(y,x); }



/// @brief Polynomial approximation of sin(x). Much faster then std libraries.
/// @tparam Order Polynomial approximation order [5,7,9]
/// @param x Angle in radians. Bound to [-pi/2,5pi/2].
/// @return Polynomial approximation of: sin(x)
template<int Order>
inline double MySin(double x)
{
    // Polynomial approximations can be found in:
    // Approximations for digital computers,
    // https://blasingame.engr.tamu.edu/z_zCourse_Archive/P620_18C/P620_zReference/PDF_Txt_Hst_Apr_Cmp_(1955).pdf

    constexpr double pi = M_PI;
    constexpr double pi_2 = M_PI_2;
    
    if      (-pi_2 <= x && x <= pi_2) {}
    else if ( pi_2 <= x && x <= 3.0 * pi_2)
        x = -(x - pi);
    else if ( 3.0 * pi_2 <= x && x <= 5.0 * pi_2)
        x = x - 2.0 * pi;
    else
        exit_on_error("sin input outside domain: " + std::to_string(x));
    
    x *= 2.0f / pi;

    if constexpr(Order == 5)
    {
        const double a1 =  1.5706268;
        const double a3 = -0.6432292;
        const double a5 =  0.0727102;
        double xx = x * x;
        return x * (a1 + xx * (a3 + xx * a5));
    }

    if constexpr(Order == 7)
    {
        const double a1 =  1.570794852;
        const double a3 = -0.645920978;
        const double a5 =  0.079487663;
        const double a7 = -0.004362476;
        double xx = x * x;
        return x * (a1 + xx * (a3 + xx * (a5 + xx * a7)));
    }

    if constexpr(Order == 9)
    {
        const double a1 =  1.57079631847;
        const double a3 = -0.64596371106;
        const double a5 =  0.07968967928;
        const double a7 = -0.00467376557;
        const double a9 =  0.00015148419;
        double xx = x * x;
        return x * (a1 + xx * (a3 + xx * (a5 + xx * (a7 + xx * a9))));
    }
}
/// @brief Polynomial order given by the macro APPROXIMATION_ORDER in ControlFlow.hh
/// @param x Angle in radians. Bound to [-pi/2,5pi/2].
/// @return Polynomial approximation of: sin(x)
inline double MySin(double x)
{ return MySin<APPROXIMATION_ORDER>(x); }



/// @brief Polynomial approximation of cos(x). Much faster then std libraries.
/// @tparam Order Polynomial approximation order [5,7,9]
/// @param x Angle in radians. Bound to [-pi/2,5pi/2].
/// @return Polynomial approximation of: cos(x)
template<int Order>
inline float MyCos(float x)
{
    constexpr float pi = M_PI;
    constexpr float pi_2 = M_PI_2;
    if (x >= pi)
        return -MySin<Order>(x - pi_2);
    else
        return MySin<Order>(x + pi_2);
}
/// @brief Polynomial order given by the macro APPROXIMATION_ORDER in ControlFlow.hh
/// @param x Angle in radians. Bound to [-pi/2,5pi/2].
/// @return Polynomial approximation of: cos(x)
inline double MyCos(double x)
{ return MyCos<APPROXIMATION_ORDER>(x); }



// Invert quaternion.
inline glm::quat Invert(const glm::quat& q)
{
    double a = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
    return glm::quat(q.w/a, -q.x/a, -q.y/a, -q.z/a);
}



// Get sign of T
template <typename T>
inline int sgn(T val)
{ return (T(0) < val) - (val < T(0)); }



// Marker for debugging
inline void Marker(std::string name="", bool newline=true)
{
    static size_t i = 0;
    std::cout << name << ": " << i << std::endl;
    if(newline)
        std::cout << std::endl;
    i++;
}



// Number of digits of given number.
inline int numDigits(int number)
{
    int digits = 0;
    while (number)
	{
        number /= 10;
        digits++;
    }
    return digits;
}



// Converts frame number to correct format, e.g. 1 -> 0001.
inline std::string FrameNumber(int f)
{
	std::string frameNumber;
	int maxDigits = 4;
	int fDigits   = numDigits(f);

	for(int i=0; i<(maxDigits-fDigits); i++)
		frameNumber += "0";
    if(f!=0)
    	frameNumber += std::to_string(f);

	return frameNumber;
}



inline std::string Format(const double n, const int precision=6)
{
    std::string output;

    // Leading space for positive numbers:
    if(n == 0 and std::signbit(n) == 0){output = " ";}
    else if(n == 0 and std::signbit(n) == 1){output = "";}
	else if(n >= 0){output = " ";}
	else   {output = "";}

    // Leading space if 2 digit number:
    //if(abs(n)<100){output += " ";}
    // Leading space if 1 digit number:
    //if(abs(n)<10){output += " ";}

    // Number of decimal digits:
    std::ostringstream ss;
    ss.precision(precision);
    ss << std::fixed << n; 
    output += ss.str();

    return output;    
}



// Print formatted double.
inline void PrintDouble(const double d, std::string name, bool newline=false, const int precision=6)
{
    std::cout << name << " = "<< Format(d,precision) << "\n";
    if(newline) std::cout << "\n";
}



// Print Quaternion:
inline void qPrint(glm::quat q, std::string name, bool newline=false, const int precision=6)
{
    std::cout << name << " = ("
    << Format(q.x,precision) << ","
    << Format(q.y,precision) << ","
    << Format(q.z,precision) << ","
    << Format(q.w,precision) << ")\n";
    if(newline) std::cout << "\n";
}
inline double qNorm(glm::quat q)
{
    return sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
}



// Print formatted Matrix.
inline void PrintMat(const Eigen::MatrixXd& A, int precision=6)
{
    for(int i=0; i<A.cols(); i++)
    {
        std::cout << "(" << Format(A(i,0),precision);
        for(int j=1; j<A.rows(); j++)
            std::cout << ", " << Format(A(i,j),precision);
        std::cout << ")" << std::endl;
    }
}



// File management:
// Creates directory in current working directory.
// Nested directories possible, e.g: path="data/output/xValues"
inline bool CreateDirectory(std::string path)
{
	namespace fs = std::filesystem;
    fs::path currentPath = fs::current_path();
	std::string directoryPath = currentPath/path;
    return fs::create_directories(directoryPath);
}
#endif //__INCLUDE_GUARD_Utility_hh__