#include "Vector2.h"



// Constructors:
Vector2::Vector2()
{ data[0] = data[1] = 0; }
Vector2::Vector2(double x, double y)
{ data[0] = x; data[1] = y; }

// Access:
double& Vector2::operator[](const int i)
{ return data[i]; }
const double& Vector2::operator[](const int i) const
{ return data[i]; }

// Basic Math:
Vector2& Vector2::operator=(const Vector2& other)
{
    data[0] = other[0];
    data[1] = other[1];
    return *this;
}
Vector2& Vector2::operator+=(const Vector2& other)
{
    data[0] += other[0];
    data[1] += other[1];
    return *this;
}
Vector2& Vector2::operator/(const double other)
{
    data[0] /= other;
    data[1] /= other;
    return *this;
}
Vector2& Vector2::operator/=(const double other)
{
    data[0] /= other;
    data[1] /= other;
    return *this;
}
Vector2 Vector2::operator+(const Vector2& other) const
{ return Vector2(data[0] + other[0], data[1] + other[1]); }
Vector2 Vector2::operator-(const Vector2& other) const
{ return Vector2(data[0] - other[0], data[1] - other[1]); }

// Boleans:
bool Vector2::operator==(const Vector2& other) const
{
    constexpr double epsilon = 1e-10;
    double diffX = data[0] - other[0];
    double diffY = data[1] - other[1];
    double sqrmag = diffX * diffX + diffY * diffY;
    return sqrmag < epsilon * epsilon;
}
bool Vector2::operator!=(const Vector2& other) const
{ return !((*this) == other); }

// Advanced Math:
double Vector2::Dot(const Vector2& a, const Vector2& b)
{ return a[0] * b[0] + a[1] * b[1]; }
double Vector2::Cross(const Vector2& a, const Vector2& b)
{ return a[0] * b[1] - b[0] * a[1]; }
double Vector2::Norm() const
{ return sqrt(Dot((*this),(*this))); }
void Vector2::Normalize()
{
    double norm = Norm();
    data[0] /= norm;
    data[1] /= norm;
}
Vector2 Vector2::Normalized() const
{
    double norm = Norm();
    return Vector2(data[0]/norm, data[1]/norm);
}

// Output:
std::ostream& operator<<(std::ostream& os, const Vector2& v) 
{
    os << v[0] << "," << v[1];
    return os;
}

// Operators that can only be defined outside of class:
Vector2 operator*(const Vector2& v, double s)
{ return Vector2(s*v[0], s*v[1]); }
Vector2 operator*(double s, const Vector2& v)
{ return Vector2(s*v[0], s*v[1]); }
Vector2 operator/(const Vector2& v, double s)
{ return Vector2(v[0]/s, v[1]/s); }