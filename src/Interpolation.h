#ifndef __INCLUDE_GUARD_Interpolation_h__
#define __INCLUDE_GUARD_Interpolation_h__
#include <math.h>           // Basic maths.
#include "TensorTypes.hh"   // General relativity tensors.
#include "Vector3.h"
#include "Vector3Int.h"
#include "Profiler.hh"      // Time measurement profiler.



// x€[0,1], f0=f(x=0), f1=f(x=1).
double LinearInterpolation(double x, double f0, double f1);
// x,y€[0,1]
// Calculates f(x,y) from f(0,0), f(0,1), f(1,0) and f(1,1):
// f(0,1) --- f(1,1)
//    |  (x,y)   |
// f(0,0) --- f(1,0)
// f10 <=> f(1,0)
double BilinearInterpolation
(double x  , double y  ,
 double f00, double f01,
 double f10, double f11);
// x,y,z[0,1]
// Calculates f(x,y,z) from f(0,0,0), f(0,0,1), f(0,1,1), ... :
double TrilinearInterpolation
(double x, double y, double z,
 double f000, double f001, double f010, double f011,
 double f100, double f101, double f110, double f111);



// x€[0,1], f0=f(x=0), f1=f(x=1), f2=f(x=2).
double SquaredInterpolation(double x, double f0, double f1, double f2);



// x€[0,1], fm1=f(-1), fp0=f(0), fp1=f(1), fp2=f(2).
double CubicInterpolation(double x, double fm1, double fp0, double fp1, double fp2);
// Cubic Interpolation, x€[0,1],
// fm1m1=f(-1,-1), fm1p0=f(-1,0), fm1p1=f(-1,1), fm1p2=f(-1,2).
// fp0m1=f( 0,-1), fp0p0=f( 0,0), fp0p1=f( 0,1), fp0p2=f( 0,2).
// fp1m1=f( 1,-1), fp1p0=f( 1,0), fp1p1=f( 1,1), fp1p2=f( 1,2).
// fp2m1=f( 2,-1), fp2p0=f( 2,0), fp2p1=f( 2,1), fp2p2=f( 2,2).
double BicubicInterpolation
(double x, double y,
 double fm1m1, double fm1p0, double fm1p1, double fm1p2,
 double fp0m1, double fp0p0, double fp0p1, double fp0p2,
 double fp1m1, double fp1p0, double fp1p1, double fp1p2,
 double fp2m1, double fp2p0, double fp2p1, double fp2p2);



bool RayTriangleIntersection(const Vector3& rayOrigin, const Vector3& rayDirection,
const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& intersectionPoint);
bool RayTriangleIntersection(const Tensor3& rayOrigin, const Tensor3& rayDirection,
const Tensor3& v0, const Tensor3& v1, const Tensor3& v2, Tensor3& intersectionPoint);



bool BarycentricWeights(const Vector3& rayOrigin, const Vector3& rayDirection,
const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& weights);
bool BarycentricWeights(const Tensor3& rayOrigin, const Tensor3& rayDirection,
const Tensor3& v0, const Tensor3& v1, const Tensor3& v2, Vector3& weights);



bool SphericalBarycentricWeights(const Tensor3& p,
const Tensor3& a, const Tensor3& b, const Tensor3& c, Tensor3& weights);
#endif //__INCLUDE_GUARD_Interpolation_h__