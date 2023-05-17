#include "Vector3.h"



// Constructors:
Vector3::Vector3()
{ data[0] = data[1] = data[2] = 0; }
Vector3::Vector3(double x, double y, double z)
{ data[0] = x; data[1] = y; data[2] = z; }

// Access:
double& Vector3::operator[](const int i)
{ return data[i]; }
const double& Vector3::operator[](const int i) const
{ return data[i]; }

// Basic Math:
Vector3& Vector3::operator=(const Vector3& other)
{
    data[0] = other[0];
    data[1] = other[1];
    data[2] = other[2];
    return *this;
}
Vector3& Vector3::operator+=(const Vector3& other)
{
    data[0] += other[0];
    data[1] += other[1];
    data[2] += other[2];
    return *this;
}
Vector3& Vector3::operator/(const double other)
{
    data[0] /= other;
    data[1] /= other;
    data[2] /= other;
    return *this;
}
Vector3& Vector3::operator/=(const double other)
{
    data[0] /= other;
    data[1] /= other;
    data[2] /= other;
    return *this;
}
Vector3 Vector3::operator+(const Vector3& other) const
{ return Vector3(data[0] + other[0], data[1] + other[1], data[2] + other[2]); }
Vector3 Vector3::operator-(const Vector3& other) const
{ return Vector3(data[0] - other[0], data[1] - other[1], data[2] - other[2]); }

// Boleans:
bool Vector3::operator==(const Vector3& other) const
{
    constexpr double epsilon = 1e-10;
    double diffX = data[0] - other[0];
    double diffY = data[1] - other[1];
    double diffZ = data[2] - other[2];
    double sqrmag = diffX * diffX + diffY * diffY + diffZ * diffZ;
    return sqrmag < epsilon * epsilon;
}
bool Vector3::operator!=(const Vector3& other) const
{ return !((*this) == other); }

// Advanced Math:
double Vector3::Dot(const Vector3& a, const Vector3& b)
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vector3 Vector3::Cross(const Vector3& a, const Vector3& b)
{
    return Vector3
    (a[1] * b[2] - b[1] * a[2],
     a[2] * b[0] - b[2] * a[0],
     a[0] * b[1] - b[0] * a[1]);
}
double Vector3::Norm() const
{ return sqrt(Dot((*this),(*this))); }
void Vector3::Normalize()
{
    double norm = Norm();
    data[0] /= norm;
    data[1] /= norm;
    data[2] /= norm;
}
Vector3 Vector3::Normalized() const
{
    double norm = Norm();
    return Vector3(data[0]/norm, data[1]/norm, data[2]/norm);
}
double Vector3::Theta() const
{ return MyAtan2(sqrt(data[0] * data[0] + data[1] * data[1]), data[2]); }
double Vector3::Phi() const
{
    double phi = MyAtan2(data[1], data[0]);
    return (phi < 0) ? phi + 2 * M_PI : phi;
}
Vector3 Vector3::GetCenter(const std::vector<Vector3>& vertices)
{
    Vector3 center(0,0,0);
    for(const Vector3& v : vertices)
        center += v;
    center /= vertices.size();
    return center;
}
bool Vector3::AreCoplanar(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3)
{
    // equation of plane: a*x + b*y + c*z + d = 0
    // analogous: n*(x-x0) = 0, where x0 is any point on the plane
    Vector3 n = Vector3::Cross(p1 - p0, p2 - p0);
    double d = -Vector3::Dot(p0, n);
    return (abs(Vector3::Dot(p3, n) + d) < 1e-8);
}

// Output:
std::ostream& operator<<(std::ostream& os, const Vector3& v) 
{
    os << v[0] << "," << v[1] << "," << v[2];
    return os;
}

// Operators that can only be defined outside of class:
Vector3 operator*(const Vector3& v, double s)
{ return Vector3(s*v[0], s*v[1], s*v[2]); }
Vector3 operator*(double s, const Vector3& v)
{ return Vector3(s*v[0], s*v[1], s*v[2]); }
Vector3 operator/(const Vector3& v, double s)
{ return Vector3(v[0]/s, v[1]/s, v[2]/s); }

// Quaternion rotation:
Vector3 operator*(const glm::quat& q, const Vector3& p)
{
    Vector3 u(q.x, q.y, q.z);
    Vector3 up  = Vector3::Cross(u, p);
    Vector3 uup = Vector3::Cross(u, up);
	return p + 2.0 * ((up * q.w) + uup);
}