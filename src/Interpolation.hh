#ifndef __INCLUDE_GUARD_Interpolation_hh__
#define __INCLUDE_GUARD_Interpolation_hh__
#include <math.h>           // Basic maths.
#include "TensorTypes.hh"   // General relativity tensors.
#include "Vector3.h"
#include "Vector3Int.h"
#include "Profiler.hh"      // Time measurement profiler.



// x€[0,1], f0=f(x=0), f1=f(x=1).
inline double LinearInterpolation(double x, double f0, double f1)
{
    return f0 * (1.0 - x) + f1 * x;
}
// x,y€[0,1]
// Calculates f(x,y) from f(0,0), f(0,1), f(1,0) and f(1,1):
// f(0,1) --- f(1,1)
//    |  (x,y)   |
// f(0,0) --- f(1,0)
// f10 <=> f(1,0)
inline double BilinearInterpolation
(double x  , double y  ,
 double f00, double f01,
 double f10, double f11)
{
    double f0 = f00 * (1.0 - x) + f10 * x;
    double f1 = f01 * (1.0 - x) + f11 * x;
    return f0 * (1.0 - y) + f1 * y;
}
// x,y,z[0,1]
// Calculates f(x,y,z) from f(0,0,0), f(0,0,1), f(0,1,1), ... :
inline double TrilinearInterpolation
(double x, double y, double z,
 double f000, double f001, double f010, double f011,
 double f100, double f101, double f110, double f111)
{
    double f00 = f000 * (1.0 - x) + f100 * x;
    double f01 = f001 * (1.0 - x) + f101 * x;
    double f10 = f010 * (1.0 - x) + f110 * x;
    double f11 = f011 * (1.0 - x) + f111 * x;
    double f0 = f00 * (1.0 - y) + f10 * y;
    double f1 = f01 * (1.0 - y) + f11 * y;
    return f0 * (1.0 - z) + f1 * z;
}



// x€[0,1], f0=f(x=0), f1=f(x=1), f2=f(x=2).
inline double SquaredInterpolation(double x, double f0, double f1, double f2)
{
    double a = (f2 + f0) / 2.0 - f0;
    double b = (f2 - f0) / 2.0;
    double c = f1;
    return (a * x + b) * x + c;
}



// x€[0,1], fm1=f(-1), fp0=f(0), fp1=f(1), fp2=f(2).
inline double CubicInterpolation(double x, double fm1, double fp0, double fp1, double fp2)
{
    double c0 = fp0;
    double c1 = 0.5 * (fp1 - fm1);
    double c2 = fm1 - 2.5 * fp0 + 2.0 * fp1 - 0.5 * fp2;
    double c3 = 0.5 * (fp2 - fm1) + 1.5 * (fp0 - fp1);
    return ((c3 * x + c2) * x + c1) * x + c0;
}
// Cubic Interpolation, x€[0,1],
// fm1m1=f(-1,-1), fm1p0=f(-1,0), fm1p1=f(-1,1), fm1p2=f(-1,2).
// fp0m1=f( 0,-1), fp0p0=f( 0,0), fp0p1=f( 0,1), fp0p2=f( 0,2).
// fp1m1=f( 1,-1), fp1p0=f( 1,0), fp1p1=f( 1,1), fp1p2=f( 1,2).
// fp2m1=f( 2,-1), fp2p0=f( 2,0), fp2p1=f( 2,1), fp2p2=f( 2,2).
inline double BicubicInterpolation
(double x, double y,
 double fm1m1, double fm1p0, double fm1p1, double fm1p2,
 double fp0m1, double fp0p0, double fp0p1, double fp0p2,
 double fp1m1, double fp1p0, double fp1p1, double fp1p2,
 double fp2m1, double fp2p0, double fp2p1, double fp2p2)
{
    double Fm1 = CubicInterpolation(x, fm1m1, fp0m1, fp1m1, fp2m1);
    double Fp0 = CubicInterpolation(x, fm1p0, fp0p0, fp1p0, fp2p0);
    double Fp1 = CubicInterpolation(x, fm1p1, fp0p1, fp1p1, fp2p1);
    double Fp2 = CubicInterpolation(x, fm1p2, fp0p2, fp1p2, fp2p2);
    return CubicInterpolation(y, Fm1, Fp0, Fp1, Fp2);
}



inline bool RayTriangleIntersection(const Vector3& rayOrigin, const Vector3& rayDirection,
const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& intersectionPoint)
{
    constexpr double epsilon = 1e-8;

    Vector3 v01 = v1 - v0;
    Vector3 v02 = v2 - v0;
    Vector3 normal = Vector3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Vector3::Dot(normal, rayDirection);
    if(abs(normalDotDirection) < epsilon)
        return false;

    // Check if triangle is behind the ray
    double d = -Vector3::Dot(normal, v0);
    double t = -(d + Vector3::Dot(normal, rayOrigin)) / normalDotDirection;
    if (t < 0)
        return false;

    // Intersection of ray with plane of triangle:
    intersectionPoint = rayOrigin + t * rayDirection;

    // Check if point is actually inside triangle:
    Vector3 c0 = intersectionPoint - v0;
    Vector3 e0 = v1 - v0;
    double V2 = Vector3::Dot(normal, Vector3::Cross(e0, c0));
    if(V2 < epsilon)
        return false;

    Vector3 c1 = intersectionPoint - v1;
    Vector3 e1 = v2 - v1;
    double V0 = Vector3::Dot(normal, Vector3::Cross(e1, c1));
    if(V0 < epsilon)
        return false;

    Vector3 c2 = intersectionPoint - v2;
    Vector3 e2 = v0 - v2;
    double V1 = Vector3::Dot(normal, Vector3::Cross(e2, c2));
    if(V1 < epsilon)
        return false;

    return true;
}
inline bool RayTriangleIntersection(const Tensor3& rayOrigin, const Tensor3& rayDirection,
const Tensor3& v0, const Tensor3& v1, const Tensor3& v2, Tensor3& intersectionPoint)
{
    constexpr double epsilon = 1e-8;

    Tensor3 v01 = v1 - v0;
    Tensor3 v02 = v2 - v0;
    Tensor3 normal = Tensor3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Tensor3::Dot(normal, rayDirection);
    if(abs(normalDotDirection) < epsilon)
        return false;

    // Check if triangle is behind the ray
    double d = -Tensor3::Dot(normal, v0);
    double t = -(d + Tensor3::Dot(normal, rayOrigin)) / normalDotDirection;
    if (t < 0)
        return false;

    // Intersection of ray with plane of triangle:
    intersectionPoint = rayOrigin + t * rayDirection;

    // Check if point is actually inside triangle:
    Tensor3 c0 = intersectionPoint - v0;
    Tensor3 e0 = v1 - v0;
    double V2 = Tensor3::Dot(normal, Tensor3::Cross(e0, c0));
    if(V2 < epsilon)
        return false;

    Tensor3 c1 = intersectionPoint - v1;
    Tensor3 e1 = v2 - v1;
    double V0 = Tensor3::Dot(normal, Tensor3::Cross(e1, c1));
    if(V0 < epsilon)
        return false;

    Tensor3 c2 = intersectionPoint - v2;
    Tensor3 e2 = v0 - v2;
    double V1 = Tensor3::Dot(normal, Tensor3::Cross(e2, c2));
    if(V1 < epsilon)
        return false;

    return true;
}



inline bool BarycentricWeights(const Vector3& rayOrigin, const Vector3& rayDirection,
const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& weights)
{
    constexpr double epsilon = 1e-8;

    Vector3 v01 = v1 - v0;
    Vector3 v02 = v2 - v0;
    Vector3 normal = Vector3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Vector3::Dot(normal, rayDirection);
    if(abs(normalDotDirection) < epsilon)
        return false;

    // Check if triangle is behind the ray
    double d = -Vector3::Dot(normal, v0);
    double t = -(d + Vector3::Dot(normal, rayOrigin)) / normalDotDirection;
    if (t < 0)
        return false;

    // Intersection of ray with plane of triangle:
    Vector3 intersectionPoint = rayOrigin + t * rayDirection;

    // Check if point is actually inside triangle:
    Vector3 c0 = intersectionPoint - v0;
    Vector3 e0 = v1 - v0;
    double V2 = Vector3::Dot(normal, Vector3::Cross(e0, c0));
    if(V2 < epsilon)
        return false;

    Vector3 c1 = intersectionPoint - v1;
    Vector3 e1 = v2 - v1;
    double V0 = Vector3::Dot(normal, Vector3::Cross(e1, c1));
    if(V0 < epsilon)
        return false;

    Vector3 c2 = intersectionPoint - v2;
    Vector3 e2 = v0 - v2;
    double V1 = Vector3::Dot(normal, Vector3::Cross(e2, c2));
    if(V1 < epsilon)
        return false;

    // Determine Barycentric Weights:
    double V = Vector3::Dot(normal, normal);
    weights[0] = V0 / V;
    weights[1] = V1 / V;
    weights[2] = V2 / V;

    return true;
}
inline bool BarycentricWeights(const Tensor3& rayOrigin, const Tensor3& rayDirection,
const Tensor3& v0, const Tensor3& v1, const Tensor3& v2, Tensor3& weights)
{
    constexpr double epsilon = 1e-8;

    Tensor3 v01 = v1 - v0;
    Tensor3 v02 = v2 - v0;
    Tensor3 normal = Tensor3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Tensor3::Dot(normal, rayDirection);
    if(abs(normalDotDirection) < epsilon)
        return false;

    // Check if triangle is behind the ray
    double d = -Tensor3::Dot(normal, v0);
    double t = -(d + Tensor3::Dot(normal, rayOrigin)) / normalDotDirection;
    if (t < 0)
        return false;

    // Intersection of ray with plane of triangle:
    Tensor3 intersectionPoint = rayOrigin + t * rayDirection;

    // Check if point is actually inside triangle:
    Tensor3 c0 = intersectionPoint - v0;
    Tensor3 e0 = v1 - v0;
    double V2 = Tensor3::Dot(normal, Tensor3::Cross(e0, c0));
    if(V2 < epsilon)
        return false;

    Tensor3 c1 = intersectionPoint - v1;
    Tensor3 e1 = v2 - v1;
    double V0 = Tensor3::Dot(normal, Tensor3::Cross(e1, c1));
    if(V0 < epsilon)
        return false;

    Tensor3 c2 = intersectionPoint - v2;
    Tensor3 e2 = v0 - v2;
    double V1 = Tensor3::Dot(normal, Tensor3::Cross(e2, c2));
    if(V1 < epsilon)
        return false;

    // Determine Barycentric Weights:
    double V = Tensor3::Dot(normal, normal);
    weights[1] = V0 / V;
    weights[2] = V1 / V;
    weights[3] = V2 / V;

    return true;
}



inline bool SphericalBarycentricWeights(const Tensor3& p,
const Tensor3& a, const Tensor3& b, const Tensor3& c, Tensor3& weights)
{
    Tensor3 vab = b - a;
    Tensor3 vbc = c - b;
    Tensor3 vca = a - c;
    
    Tensor3 vap = p - a;
    Tensor3 vbp = p - b;
    Tensor3 vcp = p - c;
    
    Tensor3 nab = Tensor3::Cross(b,a);
    Tensor3 nbc = Tensor3::Cross(c,b);
    Tensor3 nca = Tensor3::Cross(a,c);
    
    Tensor3 nap = Tensor3::Cross(p,a);
    Tensor3 nbp = Tensor3::Cross(p,b);
    Tensor3 ncp = Tensor3::Cross(p,c);

    double angleCAB = M_PI - acos(Tensor3::Dot(nca, nab));
    double angleABC = M_PI - acos(Tensor3::Dot(nab, nbc));
    double angleBCA = M_PI - acos(Tensor3::Dot(nbc, nca));

    double angleABP = M_PI - acos(Tensor3::Dot(nab, nbp));
    double angleBPA = M_PI - acos(Tensor3::Dot(nbp, -nap));
    double anglePAB = M_PI - acos(Tensor3::Dot(-nap, nab));

    double angleBCP = M_PI - acos(Tensor3::Dot(nbc, ncp));
    double angleCPB = M_PI - acos(Tensor3::Dot(ncp, -nbp));
    double anglePBC = M_PI - acos(Tensor3::Dot(-nbp, nbc));

    double angleCAP = M_PI - acos(Tensor3::Dot(nca, nap));
    double angleAPC = M_PI - acos(Tensor3::Dot(nap, -ncp));
    double anglePCA = M_PI - acos(Tensor3::Dot(-ncp, nca));

    std::cout << "angleCAB: " << angleCAB << std::endl;
    std::cout << "angleABC: " << angleABC << std::endl;
    std::cout << "angleBCA: " << angleBCA << std::endl;
    
    std::cout << "angleABP: " << angleABP << std::endl;
    std::cout << "angleBPA: " << angleBPA << std::endl;
    std::cout << "anglePAB: " << anglePAB << std::endl;
    
    std::cout << "angleBCP: " << angleBCP << std::endl;
    std::cout << "angleCPB: " << angleCPB << std::endl;
    std::cout << "anglePBC: " << anglePBC << std::endl;
    
    std::cout << "angleCAP: " << angleCAP << std::endl;
    std::cout << "angleAPC: " << angleAPC << std::endl;
    std::cout << "anglePCA: " << anglePCA << std::endl;
    
    double areaA = (angleBCP + angleCPB + anglePBC - M_PI);
    double areaB = (angleCAP + angleAPC + anglePCA - M_PI);
    double areaC = (angleABP + angleBPA + anglePAB - M_PI);
    double areaABC = areaA + areaB + areaC;
    
    
    std::cout << "areaA: " << areaA << std::endl;
    std::cout << "areaB: " << areaB << std::endl;
    std::cout << "areaC: " << areaC << std::endl;
    std::cout << "areaABC: " << areaABC << std::endl;

    return true;
}
#endif //__INCLUDE_GUARD_Interpolation_hh__