#include "Interpolation.h"

constexpr double epsilon = 1e-8;

// x€[0,1], f0=f(x=0), f1=f(x=1).
double LinearInterpolation(double x, double f0, double f1)
{
    return f0 * (1.0 - x) + f1 * x;
}
// x,y€[0,1]
// Calculates f(x,y) from f(0,0), f(0,1), f(1,0) and f(1,1):
// f(0,1) --- f(1,1)
//    |  (x,y)   |
// f(0,0) --- f(1,0)
// f10 <=> f(1,0)
double BilinearInterpolation(double x, double y,
                             double f00, double f01,
                             double f10, double f11)
{
    double f0 = f00 * (1.0 - x) + f10 * x;
    double f1 = f01 * (1.0 - x) + f11 * x;
    return f0 * (1.0 - y) + f1 * y;
}
// x,y,z[0,1]
// Calculates f(x,y,z) from f(0,0,0), f(0,0,1), f(0,1,1), ... :
double TrilinearInterpolation(double x, double y, double z,
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
double SquaredInterpolation(double x, double f0, double f1, double f2)
{
    double a = (f2 + f0) / 2.0 - f0;
    double b = (f2 - f0) / 2.0;
    double c = f1;
    return (a * x + b) * x + c;
}

// x€[0,1], fm1=f(-1), fp0=f(0), fp1=f(1), fp2=f(2).
double CubicInterpolation(double x, double fm1, double fp0, double fp1, double fp2)
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
double BicubicInterpolation(double x, double y,
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

bool RayTriangleIntersection(const Vector3 &rayOrigin, const Vector3 &rayDirection,
                             const Vector3 &v0, const Vector3 &v1, const Vector3 &v2, Vector3 &intersectionPoint)
{
    Vector3 v01 = v1 - v0;
    Vector3 v02 = v2 - v0;
    Vector3 normal = Vector3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Vector3::Dot(normal, rayDirection);
    if (abs(normalDotDirection) < epsilon)
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
    if (V2 < -epsilon)
        return false;

    Vector3 c1 = intersectionPoint - v1;
    Vector3 e1 = v2 - v1;
    double V0 = Vector3::Dot(normal, Vector3::Cross(e1, c1));
    if (V0 < -epsilon)
        return false;

    Vector3 c2 = intersectionPoint - v2;
    Vector3 e2 = v0 - v2;
    double V1 = Vector3::Dot(normal, Vector3::Cross(e2, c2));
    if (V1 < -epsilon)
        return false;

    return true;
}
bool RayTriangleIntersection(const Tensor3 &rayOrigin, const Tensor3 &rayDirection,
                             const Tensor3 &v0, const Tensor3 &v1, const Tensor3 &v2, Tensor3 &intersectionPoint)
{
    Tensor3 v01 = v1 - v0;
    Tensor3 v02 = v2 - v0;
    Tensor3 normal = Tensor3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Tensor3::Dot(normal, rayDirection);
    if (abs(normalDotDirection) < epsilon)
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
    if (V2 < -epsilon)
        return false;

    Tensor3 c1 = intersectionPoint - v1;
    Tensor3 e1 = v2 - v1;
    double V0 = Tensor3::Dot(normal, Tensor3::Cross(e1, c1));
    if (V0 < -epsilon)
        return false;

    Tensor3 c2 = intersectionPoint - v2;
    Tensor3 e2 = v0 - v2;
    double V1 = Tensor3::Dot(normal, Tensor3::Cross(e2, c2));
    if (V1 < -epsilon)
        return false;

    return true;
}

bool BarycentricWeights(const Vector3 &rayOrigin, const Vector3 &rayDirection,
                        const Vector3 &v0, const Vector3 &v1, const Vector3 &v2, Vector3 &weights)
{
    Vector3 v01 = v1 - v0;
    Vector3 v02 = v2 - v0;
    Vector3 normal = Vector3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Vector3::Dot(normal, rayDirection);
    if (abs(normalDotDirection) < epsilon)
        return false;

    // Check if triangle is behind the ray:
    double d = -Vector3::Dot(normal, v0);
    double t = -(d + Vector3::Dot(normal, rayOrigin)) / normalDotDirection;
    if (t < 0)
        return false;

    // Intersection of ray with plane of triangle:
    Vector3 p = rayOrigin + t * rayDirection;

    // Check if point is actually inside triangle:
    Vector3 c0 = p - v0;
    Vector3 e0 = v1 - v0;
    double V2 = Vector3::Dot(normal, Vector3::Cross(e0, c0));
    if (V2 < -epsilon)
        return false;

    Vector3 c1 = p - v1;
    Vector3 e1 = v2 - v1;
    double V0 = Vector3::Dot(normal, Vector3::Cross(e1, c1));
    if (V0 < -epsilon)
        return false;

    Vector3 c2 = p - v2;
    Vector3 e2 = v0 - v2;
    double V1 = Vector3::Dot(normal, Vector3::Cross(e2, c2));
    if (V1 < -epsilon)
        return false;

    // Determine Barycentric Weights:
    double V = Vector3::Dot(normal, normal);
    weights[0] = std::clamp(V0 / V, 0.0, 1.0);
    weights[1] = std::clamp(V1 / V, 0.0, 1.0);
    weights[2] = std::clamp(V2 / V, 0.0, 1.0);

    return true;
}
bool BarycentricWeights(const Tensor3 &rayOrigin, const Tensor3 &rayDirection,
                        const Tensor3 &v0, const Tensor3 &v1, const Tensor3 &v2, Vector3 &weights)
{
    Tensor3 v01 = v1 - v0;
    Tensor3 v02 = v2 - v0;
    Tensor3 normal = Tensor3::Cross(v01, v02);

    // Check if ray is orthogonal to normal:
    double normalDotDirection = Tensor3::Dot(normal, rayDirection);
    if (abs(normalDotDirection) < epsilon)
        return false;

    // Check if triangle is behind the ray:
    double d = -Tensor3::Dot(normal, v0);
    double t = -(d + Tensor3::Dot(normal, rayOrigin)) / normalDotDirection;
    if (t < 0)
        return false;

    // Intersection of ray with plane of triangle:
    Tensor3 p = rayOrigin + t * rayDirection;

    // Check if point is actually inside triangle:
    Tensor3 c0 = p - v0;
    Tensor3 e0 = v1 - v0;
    double V2 = Tensor3::Dot(normal, Tensor3::Cross(e0, c0));
    if (V2 < -epsilon)
        return false;

    Tensor3 c1 = p - v1;
    Tensor3 e1 = v2 - v1;
    double V0 = Tensor3::Dot(normal, Tensor3::Cross(e1, c1));
    if (V0 < -epsilon)
        return false;

    Tensor3 c2 = p - v2;
    Tensor3 e2 = v0 - v2;
    double V1 = Tensor3::Dot(normal, Tensor3::Cross(e2, c2));
    if (V1 < -epsilon)
        return false;

    // Determine Barycentric Weights:
    double V = Tensor3::Dot(normal, normal);
    weights[0] = std::clamp(V0 / V, 0.0, 1.0);
    weights[1] = std::clamp(V1 / V, 0.0, 1.0);
    weights[2] = std::clamp(V2 / V, 0.0, 1.0);

    return true;
}

bool SphericalBarycentricWeights(const Tensor3 &p,
                                 const Tensor3 &a, const Tensor3 &b, const Tensor3 &c, Tensor3 &weights)
{
    Tensor3 vab = b - a;
    Tensor3 vbc = c - b;
    Tensor3 vca = a - c;

    Tensor3 vap = p - a;
    Tensor3 vbp = p - b;
    Tensor3 vcp = p - c;

    Tensor3 nab = Tensor3::Cross(a, b).EuklNormalized();
    Tensor3 nbc = Tensor3::Cross(b, c).EuklNormalized();
    Tensor3 nca = Tensor3::Cross(c, a).EuklNormalized();

    Tensor3 nap = Tensor3::Cross(a, p).EuklNormalized();
    Tensor3 nbp = Tensor3::Cross(b, p).EuklNormalized();
    Tensor3 ncp = Tensor3::Cross(c, p).EuklNormalized();

    double angleCAB = M_PI - MyAcos(Tensor3::Dot(nca, nab));
    double angleABC = M_PI - MyAcos(Tensor3::Dot(nab, nbc));
    double angleBCA = M_PI - MyAcos(Tensor3::Dot(nbc, nca));

    double angleABP = M_PI - MyAcos(Tensor3::Dot(nab, nbp));
    double angleBPA = M_PI - MyAcos(Tensor3::Dot(nbp, -nap));
    double anglePAB = M_PI - MyAcos(Tensor3::Dot(-nap, nab));

    double angleBCP = M_PI - MyAcos(Tensor3::Dot(nbc, ncp));
    double angleCPB = M_PI - MyAcos(Tensor3::Dot(ncp, -nbp));
    double anglePBC = M_PI - MyAcos(Tensor3::Dot(-nbp, nbc));

    double angleCAP = M_PI - MyAcos(Tensor3::Dot(nca, nap));
    double angleAPC = M_PI - MyAcos(Tensor3::Dot(nap, -ncp));
    double anglePCA = M_PI - MyAcos(Tensor3::Dot(-ncp, nca));

    Tensor3 orientationABC = Tensor3::Cross(vab, vbc);
    Tensor3 orientationABP = Tensor3::Cross(vab, vbp);
    Tensor3 orientationBCP = Tensor3::Cross(vbc, vcp);
    Tensor3 orientationCAP = Tensor3::Cross(vca, vap);
    double areaA = (angleBCP + angleCPB + anglePBC - M_PI) * sgn(Tensor3::Dot(orientationABC, orientationBCP));
    double areaB = (angleCAP + angleAPC + anglePCA - M_PI) * sgn(Tensor3::Dot(orientationABC, orientationCAP));
    double areaC = (angleABP + angleBPA + anglePAB - M_PI) * sgn(Tensor3::Dot(orientationABC, orientationABP));
    double areaABC = areaA + areaB + areaC;

    if (areaA < -epsilon || areaB < -epsilon || areaC < -epsilon)
        return false;

    weights[1] = areaA / areaABC;
    weights[2] = areaB / areaABC;
    weights[3] = areaC / areaABC;
    return true;
}