#include "SpecialMath.h"

double Dot(const Tensor3 &vector0, const Tensor3 &vector1, const Tensor3x3 &g_ll)
{
    double dot = 0;
    for (int i = 1; i < 4; i++)
        for (int j = 1; j < 4; j++)
            dot += g_ll[{i, j}] * vector0[i] * vector1[j];
    return dot;
}
double Dot(const Tensor4 &vector0, const Tensor4 &vector1, const Tensor4x4 &gamma_ll)
{
    double dot = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            dot += gamma_ll[{i, j}] * vector0[i] * vector1[j];
    return dot;
}

double Norm2(const Tensor3 &vector, const Tensor3x3 &gamma_ll)
{
    double norm2 = 0;
    for (int i = 1; i < 4; i++)
        for (int j = 1; j < 4; j++)
            norm2 += gamma_ll[{i, j}] * vector[i] * vector[j];
    return norm2;
}

double Norm2(const Tensor4 &vector, const Tensor4x4 &g_ll)
{
    double norm2 = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            norm2 += g_ll[{i, j}] * vector[i] * vector[j];
    return norm2;
}

Tensor4 NullNormalize(const Tensor4 &vector, const Tensor4x4 &g_ll)
{
    double a = 0;
    for (int i = 1; i < 4; i++)
        for (int j = 1; j < 4; j++)
            a += g_ll[rank2Indices{i, j}] * vector[i] * vector[j];
    double b = 0;
    for (int i = 1; i < 4; i++)
        b += g_ll[rank2Indices{0, i}] * vector[0] * vector[i];
    double c = g_ll[rank2Indices{0, 0}] * vector[0] * vector[0];
    double d = -b / a + sqrt((b * b) / (a * a) - c / a);

    return Tensor4(vector[0], vector[1] * d, vector[2] * d, vector[3] * d);
}

Tensor3 TransformIFtoLF(const Tensor3 &vector, const Tensor4x4 &tetrad)
{
    return Tensor3(tetrad[{1, 1}] * vector[1] + tetrad[{1, 2}] * vector[2] + tetrad[{1, 3}] * vector[3],
                   tetrad[{2, 1}] * vector[1] + tetrad[{2, 2}] * vector[2] + tetrad[{2, 3}] * vector[3],
                   tetrad[{3, 1}] * vector[1] + tetrad[{3, 2}] * vector[2] + tetrad[{3, 3}] * vector[3]);
}
Tensor3 TransformLFtoIF(const Tensor3 &vector, const Tensor4x4 &tetradInverse)
{
    return Tensor3(tetradInverse[{1, 1}] * vector[1] + tetradInverse[{1, 2}] * vector[2] + tetradInverse[{1, 3}] * vector[3],
                   tetradInverse[{2, 1}] * vector[1] + tetradInverse[{2, 2}] * vector[2] + tetradInverse[{2, 3}] * vector[3],
                   tetradInverse[{3, 1}] * vector[1] + tetradInverse[{3, 2}] * vector[2] + tetradInverse[{3, 3}] * vector[3]);
}

Tensor4x4 TransformIFtoLF(const Tensor4x4 &tensor, const Tensor4x4 &tetrad)
{
    Tensor4x4 result(0.0);
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            for (int A = 0; A < 4; A++)
                for (int B = 0; B < 4; B++)
                    result[{a, b}] += tensor[{A, B}] * tetrad[{a, A}] * tetrad[{b, B}];
    return result;
}
Tensor4x4 TransformLFtoIF(const Tensor4x4 &tensor, const Tensor4x4 &tetradInverse)
{
    Tensor4x4 result(0.0);
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            for (int A = 0; A < 4; A++)
                for (int B = 0; B < 4; B++)
                    result[{a, b}] += tensor[{A, B}] * tetradInverse[{a, A}] * tetradInverse[{b, B}];
    return result;
}



Tensor4x4 BoostMatrix(const Tensor3 & u)
{
    double betaSq = Tensor3::Dot(u, u);
    double gamma = 1.0 / sqrt(1.0 - betaSq);
    if (betaSq < 1e-10)
        return Tensor4x4(1, 0, 0, 0, 
                         0, 1, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, 1);
    return Tensor4x4(        gamma,                         -gamma * u[1],                           -gamma * u[2],                          -gamma * u[3],
                     -gamma * u[1], 1 + (gamma - 1) * u[1] * u[1] / betaSq,     (gamma - 1) * u[1] * u[2] / betaSq,     (gamma - 1) * u[1] * u[3] / betaSq,
                     -gamma * u[2],     (gamma - 1) * u[1] * u[2] / betaSq, 1 + (gamma - 1) * u[2] * u[2] / betaSq,     (gamma - 1) * u[2] * u[3] / betaSq,
                     -gamma * u[3],     (gamma - 1) * u[1] * u[3] / betaSq,     (gamma - 1) * u[2] * u[3] / betaSq, 1 + (gamma - 1) * u[3] * u[3] / betaSq);
}



template <class FrameIn, class FrameOut>
Tensor3 Vec3ObservedByEulObs(const Tensor4 &u, const Coord &xyz, Metric &metric)
{
    if constexpr (std::is_same<FrameIn, IF>::value && std::is_same<FrameOut, IF>::value)
    { // IF -> IF
        double alpha = metric.GetAlpha(xyz);
        return Tensor3(u[1] / alpha, u[2] / alpha, u[3] / alpha);
    }
    if constexpr (std::is_same<FrameIn, IF>::value && std::is_same<FrameOut, LF>::value)
    { // IF -> LF
        double alpha = metric.GetAlpha(xyz);
        Tensor3 v(u[1] / alpha, u[2] / alpha, u[3] / alpha);
        return TransformIFtoLF(v, metric.GetTetrad(xyz));
    }
    if constexpr (std::is_same<FrameIn, LF>::value && std::is_same<FrameOut, IF>::value)
    { // LF -> IF
        double alpha = metric.GetAlpha(xyz);
        Tensor3 beta_u = metric.GetBeta_u(xyz);
        Tensor3 v((u[1] + beta_u[1]) / alpha, (u[2] + beta_u[2]) / alpha, (u[3] + beta_u[3]) / alpha);
        return TransformLFtoIF(v, metric.GetTetradInverse(xyz));
    }
    if constexpr (std::is_same<FrameIn, LF>::value && std::is_same<FrameOut, LF>::value)
    { // LF -> LF
        double alpha = metric.GetAlpha(xyz);
        Tensor3 beta_u = metric.GetBeta_u(xyz);
        return Tensor3((u[1] + beta_u[1]) / alpha, (u[2] + beta_u[2]) / alpha, (u[3] + beta_u[3]) / alpha);
    }
}

template Tensor3 Vec3ObservedByEulObs<IF, IF>(const Tensor4 &u, const Coord &xyz, Metric &metric);
template Tensor3 Vec3ObservedByEulObs<IF, LF>(const Tensor4 &u, const Coord &xyz, Metric &metric);
template Tensor3 Vec3ObservedByEulObs<LF, IF>(const Tensor4 &u, const Coord &xyz, Metric &metric);
template Tensor3 Vec3ObservedByEulObs<LF, LF>(const Tensor4 &u, const Coord &xyz, Metric &metric);

Vector3 CircumCenter(const Vector3 &a, const Vector3 &b, const Vector3 &c)
{
    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 normal = Vector3::Cross(ab, ac).Normalized();
    Vector3 e1 = ab.Normalized();            // x-axes of 2d plane in 3d space.
    Vector3 e2 = Vector3::Cross(normal, e1); // y-axes of 2d plane in 3d space.
    // Planar projection of a,b,c onto plane of triangle:
    Vector2 a2 = Vector2(0, 0);
    Vector2 b2 = Vector2(Vector3::Dot(e1, ab), Vector3::Dot(e2, ab));
    Vector2 c2 = Vector2(Vector3::Dot(e1, ac), Vector3::Dot(e2, ac));
    // Find center and convert back to 3d space:
    Vector2 center = CircumCenter(a2, b2, c2);
    return center[0] * e1 + center[1] * e2 + a;
}
Vector2 CircumCenter(const Vector2 &a, const Vector2 &b, const Vector2 &c)
{
    Vector2 dirAB = (b - a) / 2.0;
    Vector2 dirAC = (c - a) / 2.0;
    Vector2 midAB = a + dirAB;
    Vector2 midAC = a + dirAC;
    Vector2 perpAB = Vector2(-dirAB[1], dirAB[0]);
    Vector2 perpAC = Vector2(-dirAC[1], dirAC[0]);
    return LineLineIntersection(midAB, perpAB, midAC, perpAC);
}

double PointToPointDistance(const Vector3 &point0, const Vector3 &point1)
{
    return (point1 - point0).Norm();
}
bool IsParallel(const Vector3 &direction0, const Vector3 &direction1)
{
    return 1.0 - abs(Vector3::Dot(direction0.Normalized(), direction1.Normalized())) < 1e-4;
}
Vector2 LineLineIntersection(const Vector2 &lineSupport0, const Vector2 &lineDirection0, const Vector2 &lineSupport1, const Vector2 &lineDirection1)
{
    double t = Vector2::Cross(lineSupport1 - lineSupport0, lineDirection1) / Vector2::Cross(lineDirection0, lineDirection1);
    return lineSupport0 + (t * lineDirection0);
}
Vector3 LinePlaneIntersection(const Vector3 &lineSupport, const Vector3 &lineDirection, const Vector3 &planeSupport, const Vector3 &planeNormal)
{
    double d = Vector3::Dot(planeSupport, planeNormal);
    double t = (d - Vector3::Dot(lineSupport, planeNormal)) / Vector3::Dot(lineDirection, planeNormal);
    return lineSupport + t * lineDirection;
}
Vector3 LinePlaneIntersection(const Vector3 &lineSupport, const Vector3 &lineDirection, const Vector3 &planeSupport, const Vector3 &planeNormal, double &t)
{
    double d = Vector3::Dot(planeSupport, planeNormal);
    t = (d - Vector3::Dot(lineSupport, planeNormal)) / Vector3::Dot(lineDirection, planeNormal);
    return lineSupport + t * lineDirection;
}

std::vector<Vector3> SortPointsOnSphereAroundCenter(const std::vector<Vector3> &points)
{
    Vector3 center = Vector3(0, 0, 0);
    for (const Vector3 &p : points)
        center += p;
    center = center.Normalized();

    // Rotate center to (th,ph) = (pi/2,pi) and turn to spherical coordinates:
    glm::vec3 from(center[0], center[1], center[2]);
    glm::vec3 to(-1, 0, 0);
    glm::quat q = glm::quat(from, to);
    double thetaCenter = M_PI / 2.0;
    double phiCenter = M_PI;

    std::vector<double> angles;
    angles.reserve(points.size());
    for (const Vector3 &p : points)
    {
        // Angular position in localized spherical coordinates:
        Vector3 rotated = q * p;
        double theta = rotated.Theta();
        double phi = rotated.Phi();

        // Angles of points anti-clockwise around center vertex:
        double angle = MyAtan2(theta - thetaCenter, phi - phiCenter);
        angles.push_back(angle);
    }

    // Sort points to be anti-clockwise around its center:
    std::vector<int> indexArray;
    indexArray.reserve(points.size());
    for (int i = 0; i < points.size(); i++)
        indexArray.push_back(i);

    std::sort(indexArray.begin(), indexArray.end(),
              [&angles](int i, int j)
              { return angles[i] < angles[j]; });

    std::vector<Vector3> sortedPoints;
    sortedPoints.reserve(points.size());
    for (int i = 0; i < points.size(); i++)
        sortedPoints.push_back(points[indexArray[i]]);

    return sortedPoints;
}

double SphericalTriangleArea(const Vector3 &a, const Vector3 &b, const Vector3 &c, double radius)
{ // https://www.johndcook.com/blog/2021/11/29/area-of-spherical-triangle/
    double A = abs(Vector3::Dot(a, Vector3::Cross(b, c)));
    double B = 1.0 + Vector3::Dot(a, b) + Vector3::Dot(b, c) + Vector3::Dot(c, a);
    double E = 2.0 * MyAtan(A / B);
    return radius * radius * E;
}