#ifndef __INCLUDE_GUARD_SpecialMath_h__
#define __INCLUDE_GUARD_SpecialMath_h__
#include <vector>                     // std::vector<T>
#include <algorithm>                  // sort
#include "glm/glm/gtc/quaternion.hpp" // Quaternions.
#include "DataTypes.hh"               // General relativity tensors.
#include "Spacetimes.h"               // Metric data.
#include "Vector2.h"                  // double Vector2 type.
#include "Vector3.h"                  // double Vector3 type.

double Dot(const Tensor3 &vector0, const Tensor3 &vector1, const Tensor3x3 &g_ll);
double Dot(const Tensor4 &vector0, const Tensor4 &vector1, const Tensor4x4 &gamma_ll);
double Norm2(const Tensor3 &vector, const Tensor3x3 &gamma_ll);
double Norm2(const Tensor4 &vector, const Tensor4x4 &g_ll);
Tensor4 NullNormalize(const Tensor4 &vector, const Tensor4x4 &g_ll);

Tensor3 TransformIFtoLF(const Tensor3 &vector, const Tensor4x4 &tetrad);
Tensor3 TransformLFtoIF(const Tensor3 &vector, const Tensor4x4 &tetradInverse);

Tensor4x4 TransformIFtoLF(const Tensor4x4 &tensor, const Tensor4x4 &tetrad);
Tensor4x4 TransformLFtoIF(const Tensor4x4 &tensor, const Tensor4x4 &tetradInverse);

template <class FrameIn, class FrameOut>
Tensor3 Vec3ObservedByEulObs(const Tensor4 &u, const Coord &xyz, Metric &metric);

Vector3 CircumCenter(const Vector3 &a, const Vector3 &b, const Vector3 &c);
Vector2 CircumCenter(const Vector2 &a, const Vector2 &b, const Vector2 &c);
Vector2 LineIntersection(const Vector2 &p1, const Vector2 &d1, const Vector2 &p2, const Vector2 &d2);

bool RayPlaneIntersection(const Vector3 &p0, const Vector3 &d, const Vector3 &p1, const Vector3 &n, Vector3 &pIntersect, double &t);

std::vector<Vector3> SortPointsOnSphereAroundCenter(const std::vector<Vector3> &points);

double SphericalTriangleArea(const Vector3 &a, const Vector3 &b, const Vector3 &c, double radius = 1.0);
#endif //__INCLUDE_GUARD_SpecialMath_h__