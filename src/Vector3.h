#ifndef __INCLUDE_GUARD_Vector3_h__
#define __INCLUDE_GUARD_Vector3_h__
#include <iostream>
#include <vector>
#include <math.h>
#include "Utility.hh"
#include "glm/glm/gtc/quaternion.hpp"

struct Vector3
{
private:
    double data[3];

public:
    // Constructors:
    Vector3();
    Vector3(double x, double y, double z);

    // Access:
    double &operator[](const int i);
    const double &operator[](const int i) const;

    // Basic Math:
    Vector3 &operator=(const Vector3 &other);
    Vector3 &operator+=(const Vector3 &other);
    Vector3 &operator/(const double other);
    Vector3 &operator/=(const double other);
    Vector3 operator+(const Vector3 &other) const;
    Vector3 operator-(const Vector3 &other) const;

    // Boleans:
    bool operator==(const Vector3 &other) const;
    bool operator!=(const Vector3 &other) const;

    // Advanced Math:
    static double Dot(const Vector3 &a, const Vector3 &b);
    static Vector3 Cross(const Vector3 &a, const Vector3 &b);
    double Norm() const;
    void Normalize();
    Vector3 Normalized() const;
    double Theta() const;
    double Phi() const;
    static Vector3 GetCenter(const std::vector<Vector3> &vertices);
    static bool AreCoplanar(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3);

    // Output:
    friend std::ostream &operator<<(std::ostream &os, const Vector3 &v);
};
std::ostream &operator<<(std::ostream &os, const Vector3 &v);

// Operators that can only be defined outside of class:
Vector3 operator*(const Vector3 &v, double s);
Vector3 operator*(double s, const Vector3 &v);
Vector3 operator/(const Vector3 &v, double s);

// Quaternion rotation:
Vector3 operator*(const glm::quat &q, const Vector3 &p);
#endif //__INCLUDE_GUARD_Vector3_h__