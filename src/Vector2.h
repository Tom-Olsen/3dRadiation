#ifndef __INCLUDE_GUARD_Vector2_h__
#define __INCLUDE_GUARD_Vector2_h__
#include <iostream>
#include <vector>
#include <math.h>



struct Vector2
{
private:
    double data[2];
public:
    // Constructors:
    Vector2();
    Vector2(double x, double y);

    // Access:
    double& operator[](const int i);
    const double& operator[](const int i) const;

    // Basic Math:
    Vector2& operator=(const Vector2& other);
    Vector2& operator+=(const Vector2& other);
    Vector2& operator/(const double other);
    Vector2& operator/=(const double other);
    Vector2 operator+(const Vector2& other) const;
    Vector2 operator-(const Vector2& other) const;
    
    // Boleans:
    bool operator==(const Vector2& other) const;
    bool operator!=(const Vector2& other) const;

    // Advanced Math:
    static double Dot(const Vector2& a, const Vector2& b);
    static double Cross(const Vector2& a, const Vector2& b);
    double Norm() const; 
    void Normalize();
    Vector2 Normalized() const;
    
    // Output:
    friend std::ostream& operator<<(std::ostream& os, const Vector2& v);
};
std::ostream& operator<<(std::ostream& os, const Vector2& v);

// Operators that can only be defined outside of class:
Vector2 operator*(const Vector2& v, double s);
Vector2 operator*(double s, const Vector2& v);
Vector2 operator/(const Vector2& v, double s);
#endif //__INCLUDE_GUARD_Vector2_h__