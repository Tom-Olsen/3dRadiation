#ifndef __INCLUDE_GUARD_TesorTypes_hh__
#define __INCLUDE_GUARD_TesorTypes_hh__
#include <iostream>                     // Output to terminal.
#include "glm/glm/gtc/quaternion.hpp"   // Quaternions.
#include "ControlFlow.hh"               // Template arguments and profiling macros.
#include "Utility.hh"                   // Utility functions.
#include "eigen/Eigen/Dense"            // Eigen library for solving linear systems.



/*
 * Access pattern for all tensors follows the mathematical access pattern:
 * A[{i,j}] = A_ij
 * 
 * This allows for mathematical equations to be translated one to one:
 * C = A*B <=> C_ij = A_ik B_kj         sum over k
 * => C[{i,j}] = A[{i,k}] * B[{k,j}]    sum over k
 * 
 * This leads to the last index being the fastest, => for(i...)for(j...)for(k...)
 *     (B00,B01,B02,B03) 
 * B = (B10,B11,B12,B13) <=> Memory of B = (B00, B01, B02, B10, B11, B12, B20, B21, B22)
 *     (B20,B21,B22,B23) 
 *     (B30,B31,B32,B33) 
 * 
 * Initialization is row major:
 *     ( 0, 1, 2, 3) 
 * B = ( 4, 5, 6, 7) <=> Tensor3x3 B(0,1,2,3, 4,5,6,7, 8,9,10,11, 12,13,14,15);
 *     ( 8, 9,10,11) 
 *     (12,13,14,15) 
 * 
 * Tensor3 has indexes 1,2,3    <=> Tensor3 a(1,2,3), a[1]=1, a[2]=2, a[3]=3, a[0] => invalid
 * Tensor4 has indexes 0,1,2,3  <=> Tensor4 a(1,2,3,4), a[0]=1, a[1]=2, a[2]=3, a[3]=4
 * 
 * Tensor3x3 and Tensor4x4 analogous, but in 2d!
 */



struct rank2Indices
{int i, j;};
struct rank3Indices
{int i, j, k;};



struct Coord
{
private:
    double data[3];
public:
    // Constructors:
    Coord(double value=0)
    { data[0] = data[1] = data[2] = value; }
    Coord(double x, double y, double z)
    { data[0] = x; data[1] = y; data[2] = z; }

    // Accessors:
    double& operator[](const int index)
    { return data[index - 1]; }
    const double& operator[](const int index) const
    { return data[index - 1]; }

    // Lengths and Angles:
    double EuklNormSquared() const
    { return data[0] * data[0] + data[1] * data[1] + data[2] * data[2]; }
    double EuklNorm() const
    { return sqrt(EuklNormSquared()); }
    double Theta()
    { return MyAtan2(sqrt(data[0] * data[0] + data[1] * data[1]), data[2]); }
    double Phi()
    {
        double phi = MyAtan2(data[1], data[0]);
        return (phi < 0) ? phi + 2 * M_PI : phi;
    }

    // Basic Math:
    static double Dot(const Coord& a, const Coord& b)
    { return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }
    static Coord Cross(const Coord& a, const Coord& b)
    {
        return Coord
        (a[2] * b[3] - b[2] * a[3],
         a[3] * b[1] - b[3] * a[1],
         a[1] * b[2] - b[1] * a[2]);
    }

    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0], precision) << ","
    	<< Format(data[1], precision) << ","
    	<< Format(data[2], precision) << ")\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Coord& x);
};
// Addition/Substraction:
inline Coord operator+(const Coord& lhs, const Coord& rhs)
{ return Coord(lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]); }
inline Coord operator-(const Coord& lhs, const Coord& rhs)
{ return Coord(lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]); }

// Scalar multiplication:
inline Coord operator*(double lhs, const Coord& rhs)
{ return Coord(lhs * rhs[1], lhs * rhs[2], lhs * rhs[3]); }
inline Coord operator*(const Coord& lhs, double rhs)
{ return Coord(lhs[1] * rhs, lhs[2] * rhs, lhs[3] * rhs); }

// Quaternion multiplication:
inline Coord operator*(const glm::quat& q, const Coord& p)
{
    Coord u(q.x, q.y, q.z);
    Coord up  = Coord::Cross(u, p);
    Coord uup = Coord::Cross(u, up);
	return p + 2.0 * ((up * q.w) + uup);
}

// Output:
inline std::ostream& operator<<(std::ostream& os, const Coord& x) 
{
    os << x[1] << "," << x[2] << "," << x[3];
    return os;
}





struct Tensor2
{
private:
    double data[2];
public:
    // Constructors:
    Tensor2(double value=0)
    { data[0] = data[1]; }
    Tensor2(double y, double z)
    { data[0] = y; data[1] = z; }
    
    // Accessors:
    double& operator[](const int index)
    { return data[index - 2]; }
    const double& operator[](const int index) const
    { return data[index - 2]; }
    
    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0], precision) << ","
    	<< Format(data[1], precision) << ")\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor2& x);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor2& x) 
{
    os << x[2] << "," << x[3];
    return os;
}





struct Tensor3
{
private:
    double data[3];
public:
    // Constructors:
    Tensor3(double value=0)
    { data[0] = data[1] = data[2] = value; }
    Tensor3(double x, double y, double z)
    { data[0] = x; data[1] = y; data[2] = z; }
    
    // Accessors:
    double& operator[](const int index)
    { return data[index - 1]; }
    const double& operator[](const int index) const
    { return data[index - 1]; }

    // Lengths and Angles:
    double EuklNorm() const
    { return sqrt(Dot((*this), (*this))); }
    Tensor3 EuklNormalized() const
    {
        double norm = EuklNorm();
        return Tensor3(data[0] / norm, data[1] / norm, data[2] / norm);
    }
    // Arc distance between a and b on unit sphere. (assumes a and b are normalized)
    static double UnitSphereNorm(const Tensor3& a, const Tensor3& b)
    {
        double crossNorm = Tensor3::Cross(a, b).EuklNorm();
        double dot = Tensor3::Dot(a, b);
        return MyAtan2(crossNorm,dot);
    }
    double Theta() const
    { return MyAtan2(sqrt(data[0] * data[0] + data[1] * data[1]), data[2]); }
    double Phi() const
    {
        double phi = MyAtan2(data[1], data[0]);
        return (phi < 0) ? phi + 2 * M_PI : phi;
    }

    // Basic Math:
    static double Dot(const Tensor3& a, const Tensor3& b)
    { return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }
    static Tensor3 Cross(const Tensor3& a, const Tensor3& b)
    {
        return Tensor3
        (a[2] * b[3] - b[2] * a[3],
         a[3] * b[1] - b[3] * a[1],
         a[1] * b[2] - b[1] * a[2]);
    }
    
    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0], precision) << ","
    	<< Format(data[1], precision) << ","
    	<< Format(data[2], precision) << ")\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor3& v);
};
// Addition/Substraction:
inline Tensor3 operator+(const Tensor3& lhs, const Tensor3& rhs)
{ return Tensor3(lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]); }
inline Tensor3 operator-(const Tensor3& lhs, const Tensor3& rhs)
{ return Tensor3(lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]); }
inline Tensor3 operator-(const Tensor3& tensor)
{ return {-tensor[1], -tensor[2], -tensor[3]}; }

// Scalar multiplication:
inline Tensor3 operator*(double lhs, const Tensor3& rhs)
{ return Tensor3(lhs * rhs[1], lhs * rhs[2], lhs * rhs[3]); }
inline Tensor3 operator*(const Tensor3& lhs, double rhs)
{ return Tensor3(lhs[1] * rhs, lhs[2] * rhs, lhs[3] * rhs); }

// Quaternion rotation:
inline Tensor3 operator*(const glm::quat& q, const Tensor3& p)
{
    Tensor3 u(q.x, q.y, q.z);
    Tensor3 up  = Tensor3::Cross(u, p);
    Tensor3 uup = Tensor3::Cross(u, up);
	return p + 2.0 * ((up * q.w) + uup);
}

// Output
inline std::ostream& operator<<(std::ostream& os, const Tensor3& v) 
{
    os << v[1] << "," << v[2] << "," << v[3];
    return os;
}





struct Tensor4
{
private:
    double data[4];
public:
    // Constructors:
    Tensor4(double value=0)
    { data[0] = data[1] = data[2] = data[3] = value; }
    Tensor4(double t, double x, double y, double z)
    { data[0] = t; data[1] = x; data[2] = y; data[3] = z; }

    // Accessors:
    double& operator[](const int index)
    { return data[index]; }
    const double& operator[](const int index) const
    { return data[index]; }
    
    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0], precision) << ","
    	<< Format(data[1], precision) << ","
    	<< Format(data[2], precision) << ","
    	<< Format(data[3], precision) << ")\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor4& v);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor4& v) 
{
    os << v[0] << "," << v[1] << "," << v[2] << "," << v[3];
    return os;
}





struct Tensor2x2
{
private:
    double data[4];
public:
    // Constructors:
    Tensor2x2(double value=0)
    { data[0] = data[1] = data[2] = data[3] = value; }
    Tensor2x2(double yy, double yz,
              double zy, double zz)
    {
        data[0*2 + 0] = yy; data[0*2 + 1] = yz;
        data[1*2 + 0] = zy; data[1*2 + 1] = zz;
    }

    // Accessors:
    double& operator[](const rank2Indices& index)
    { return data[(index.i - 2) * 2 + (index.j - 2)]; }
    const double& operator[](const rank2Indices& index) const
    { return data[(index.i - 2) * 2 + (index.j - 2)]; }
    
    // Basic Math:
    Tensor2x2 Invert()
    {
        Tensor2x2 invers;
        using namespace Eigen;
        Map<Matrix<double, 2, 2, RowMajor>> matrix(data);
        Map<Matrix<double, 2, 2, RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }

    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size,' ');
        std::cout << name  << " = (" << Format(data[0 * 2 + 0], precision) << "," << Format(data[0 * 1 + 1], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[1 * 2 + 0], precision) << "," << Format(data[1 * 1 + 1], precision) << ")" << "\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor2x2& v);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor2x2& v) 
{
    os << "(" << v[{2,2}] << "," << v[{2,3}] << "|" << v[{3,2}] << "," << v[{3,3}] << ")";
    return os;
}





struct Tensor3x3
{
private:
    double data[9];
public:
    // Constructors:
    Tensor3x3(double value=0)
    {
        for(int i=0; i<9; i++)
            data[i] = value;
    }
    Tensor3x3(double xx, double xy, double xz,
              double yx, double yy, double yz,
              double zx, double zy, double zz)
    {
        data[0] = xx; data[1] = xy; data[2] = xz;
        data[3] = yx; data[4] = yy; data[5] = yz;
        data[6] = zx; data[7] = zy; data[8] = zz;
    }

    // Accessors:
    double& operator[](const rank2Indices& index)
    { return data[(index.i - 1) * 3 + (index.j - 1)]; }
    const double& operator[](const rank2Indices& index) const
    { return data[(index.i - 1) * 3 + (index.j - 1)]; }

    // Basic Math:
    Tensor3x3 Invert()
    {
        Tensor3x3 invers;
        using namespace Eigen;
        Map<Matrix<double, 3, 3, RowMajor>> matrix(data);
        Map<Matrix<double, 3, 3, RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }
    
    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size,' ');
        std::cout << space << "   (" << Format(data[0 * 3 + 0], precision) << "," << Format(data[0 * 3 + 1], precision) << "," << Format(data[0 * 3 + 2], precision) << ")" << "\n";
        std::cout << name  << " = (" << Format(data[1 * 3 + 0], precision) << "," << Format(data[1 * 3 + 1], precision) << "," << Format(data[1 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 3 + 0], precision) << "," << Format(data[2 * 3 + 1], precision) << "," << Format(data[2 * 3 + 2], precision) << ")" << "\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor3x3& v);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor3x3& v) 
{
    os << "(" << v[{1,1}] << "," << v[{1,2}] << "," << v[{1,3}]
       << "|" << v[{2,1}] << "," << v[{2,2}] << "," << v[{2,3}]
       << "|" << v[{3,1}] << "," << v[{3,2}] << "," << v[{3,3}] << ")";
    return os;
}





struct Tensor4x4
{
private:
    double data[16];
public:
    // Constructors:
    Tensor4x4(double value=0)
    {
        for(int i=0; i<16; i++)
            data[i] = value;
    }
    Tensor4x4(double tt, double tx, double ty, double tz,
              double xt, double xx, double xy, double xz,
              double yt, double yx, double yy, double yz,
              double zt, double zx, double zy, double zz)
    {
        data[ 0] = tt; data[ 1] = tx; data[ 2] = ty; data[ 3] = tz;
        data[ 4] = xt; data[ 5] = xx; data[ 6] = xy; data[ 7] = xz;
        data[ 8] = yt; data[ 9] = yx; data[10] = yy; data[11] = yz;
        data[12] = zt; data[13] = zx; data[14] = zy; data[15] = zz;
    }

    // Accessors:
    double& operator[](const rank2Indices& index)
    { return data[index.i * 4 + index.j]; }
    const double& operator[](const rank2Indices& index) const
    { return data[index.i * 4 + index.j]; }

    // Basic Math:
    Tensor4x4 Invert()
    {
        Tensor4x4 invers;
        using namespace Eigen;
        Map<Matrix<double, 4, 4, RowMajor>> matrix(data);
        Map<Matrix<double, 4, 4, RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }
    
    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size,' ');
        std::cout << space << "   (" << Format(data[0 * 4 + 0], precision) << "," << Format(data[0 * 4 + 1], precision) << "," << Format(data[0 * 4 + 2], precision) << "," << Format(data[0 * 4 + 3], precision) << ")" << "\n";
        std::cout << name  << " = (" << Format(data[1 * 4 + 0], precision) << "," << Format(data[1 * 4 + 1], precision) << "," << Format(data[1 * 4 + 2], precision) << "," << Format(data[1 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 4 + 0], precision) << "," << Format(data[2 * 4 + 1], precision) << "," << Format(data[2 * 4 + 2], precision) << "," << Format(data[2 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[3 * 4 + 0], precision) << "," << Format(data[3 * 4 + 1], precision) << "," << Format(data[3 * 4 + 2], precision) << "," << Format(data[3 * 4 + 3], precision) << ")" << "\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor4x4& v);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor4x4& v) 
{
    os << "(" << v[{0,0}] << "," << v[{0,1}] << "," << v[{0,2}] << "," << v[{0,3}]
       << "|" << v[{1,0}] << "," << v[{1,1}] << "," << v[{1,2}] << "," << v[{1,3}]
       << "|" << v[{2,0}] << "," << v[{2,1}] << "," << v[{2,2}] << "," << v[{2,3}]
       << "|" << v[{3,0}] << "," << v[{3,1}] << "," << v[{3,2}] << "," << v[{3,3}] << ")";
    return os;
}





struct Tensor3x3x3
{
private:
    double data[27];
public:
    // Constructors:
    Tensor3x3x3(double value=0)
    {
        for(int i=0; i<27; i++)
            data[i] = value;
    }
    Tensor3x3x3
    (double xxx, double xxy, double xxz,
     double xyx, double xyy, double xyz,
     double xzx, double xzy, double xzz,
        double yxx, double yxy, double yxz,
        double yyx, double yyy, double yyz,
        double yzx, double yzy, double yzz,
            double zxx, double zxy, double zxz,
            double zyx, double zyy, double zyz,
            double zzx, double zzy, double zzz)
    {
        data[0 * 9 + 0 * 3 + 0] = xxx; data[0 * 9 + 0 * 3 + 1] = xxy; data[0 * 9 + 0 * 3 + 2] = xxz;
        data[0 * 9 + 1 * 3 + 0] = xyx; data[0 * 9 + 1 * 3 + 1] = xyy; data[0 * 9 + 1 * 3 + 2] = xyz;
        data[0 * 9 + 2 * 3 + 0] = xzx; data[0 * 9 + 2 * 3 + 1] = xzy; data[0 * 9 + 2 * 3 + 2] = xzz;
        data[1 * 9 + 0 * 3 + 0] = yxx; data[1 * 9 + 0 * 3 + 1] = yxy; data[1 * 9 + 0 * 3 + 2] = yxz;
        data[1 * 9 + 1 * 3 + 0] = yyx; data[1 * 9 + 1 * 3 + 1] = yyy; data[1 * 9 + 1 * 3 + 2] = yyz;
        data[1 * 9 + 2 * 3 + 0] = yzx; data[1 * 9 + 2 * 3 + 1] = yzy; data[1 * 9 + 2 * 3 + 2] = yzz;
        data[2 * 9 + 0 * 3 + 0] = zxx; data[2 * 9 + 0 * 3 + 1] = zxy; data[2 * 9 + 0 * 3 + 2] = zxz;
        data[2 * 9 + 1 * 3 + 0] = zyx; data[2 * 9 + 1 * 3 + 1] = zyy; data[2 * 9 + 1 * 3 + 2] = zyz;
        data[2 * 9 + 2 * 3 + 0] = zzx; data[2 * 9 + 2 * 3 + 1] = zzy; data[2 * 9 + 2 * 3 + 2] = zzz;
    }
    
    // Accessors:
    double& operator[](const rank3Indices& index)
    { return data[(index.i - 1) * 9 + (index.j - 1) * 3 + (index.k - 1)]; }
    const double& operator[](const rank3Indices& index) const
    { return data[(index.i - 1) * 9 + (index.j - 1) * 3 + (index.k - 1)]; }
    
    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size,' ');
        std::cout << space << "   (" << Format(data[0 * 9 + 0 * 3 + 0], precision) << "," << Format(data[0 * 9 + 0 * 3 + 1], precision) << "," << Format(data[0 * 9 + 0 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[0 * 9 + 1 * 3 + 0], precision) << "," << Format(data[0 * 9 + 1 * 3 + 1], precision) << "," << Format(data[0 * 9 + 1 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[0 * 9 + 2 * 3 + 0], precision) << "," << Format(data[0 * 9 + 2 * 3 + 1], precision) << "," << Format(data[0 * 9 + 2 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (--------------------------------)\n";
        std::cout << space << "   (" << Format(data[1 * 9 + 0 * 3 + 0], precision) << "," << Format(data[1 * 9 + 0 * 3 + 1], precision) << "," << Format(data[1 * 9 + 0 * 3 + 2], precision) << ")" << "\n";
        std::cout << name  << " = (" << Format(data[1 * 9 + 1 * 3 + 0], precision) << "," << Format(data[1 * 9 + 1 * 3 + 1], precision) << "," << Format(data[1 * 9 + 1 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[1 * 9 + 2 * 3 + 0], precision) << "," << Format(data[1 * 9 + 2 * 3 + 1], precision) << "," << Format(data[1 * 9 + 2 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (--------------------------------)\n";
        std::cout << space << "   (" << Format(data[2 * 9 + 0 * 3 + 0], precision) << "," << Format(data[2 * 9 + 0 * 3 + 1], precision) << "," << Format(data[2 * 9 + 0 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 9 + 1 * 3 + 0], precision) << "," << Format(data[2 * 9 + 1 * 3 + 1], precision) << "," << Format(data[2 * 9 + 1 * 3 + 2], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 9 + 2 * 3 + 0], precision) << "," << Format(data[2 * 9 + 2 * 3 + 1], precision) << "," << Format(data[2 * 9 + 2 * 3 + 2], precision) << ")" << "\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor3x3x3& v);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor3x3x3& v) 
{
    os << "(" << v[{1,1,1}] << "," << v[{1,1,2}] << "," << v[{1,1,3}]
       << "|" << v[{1,2,1}] << "," << v[{1,2,2}] << "," << v[{1,2,3}]
       << "|" << v[{1,3,1}] << "," << v[{1,3,2}] << "," << v[{1,3,3}] << ")\n";
    os << "(" << v[{2,1,1}] << "," << v[{2,1,2}] << "," << v[{2,1,3}]
       << "|" << v[{2,2,1}] << "," << v[{2,2,2}] << "," << v[{2,2,3}]
       << "|" << v[{2,3,1}] << "," << v[{2,3,2}] << "," << v[{2,3,3}] << ")\n";
    os << "(" << v[{3,1,1}] << "," << v[{3,1,2}] << "," << v[{3,1,3}]
       << "|" << v[{3,2,1}] << "," << v[{3,2,2}] << "," << v[{3,2,3}]
       << "|" << v[{3,3,1}] << "," << v[{3,3,2}] << "," << v[{3,3,3}] << ")";
    return os;
}





struct Tensor4x4x4
{
private:
    double data[64];
public:
    // Constructors:
    Tensor4x4x4(double value=0)
    {
        for(int i=0; i<64; i++)
            data[i] = value;
    }
    Tensor4x4x4
    (double ttt, double ttx, double tty, double ttz,
     double txt, double txx, double txy, double txz,
     double tyt, double tyx, double tyy, double tyz,
     double tzt, double tzx, double tzy, double tzz,
     double xtt, double xtx, double xty, double xtz,
     double xxt, double xxx, double xxy, double xxz,
     double xyt, double xyx, double xyy, double xyz,
     double xzt, double xzx, double xzy, double xzz,
     double ytt, double ytx, double yty, double ytz,
     double yxt, double yxx, double yxy, double yxz,
     double yyt, double yyx, double yyy, double yyz,
     double yzt, double yzx, double yzy, double yzz,
     double ztt, double ztx, double zty, double ztz,
     double zxt, double zxx, double zxy, double zxz,
     double zyt, double zyx, double zyy, double zyz,
     double zzt, double zzx, double zzy, double zzz)
    {
        data[0 * 16 + 0 * 4 + 0] = ttt; data[0 * 16 + 0 * 4 + 1] = ttx; data[0 * 16 + 0 * 4 + 2] = tty; data[0 * 16 + 0 * 4 + 3] = ttz;
        data[0 * 16 + 1 * 4 + 0] = txt; data[0 * 16 + 1 * 4 + 1] = txx; data[0 * 16 + 1 * 4 + 2] = txy; data[0 * 16 + 1 * 4 + 3] = txz;
        data[0 * 16 + 2 * 4 + 0] = tyt; data[0 * 16 + 2 * 4 + 1] = tyx; data[0 * 16 + 2 * 4 + 2] = tyy; data[0 * 16 + 2 * 4 + 3] = tyz;
        data[0 * 16 + 3 * 4 + 0] = tzt; data[0 * 16 + 3 * 4 + 1] = tzx; data[0 * 16 + 3 * 4 + 2] = tzy; data[0 * 16 + 3 * 4 + 3] = tzz;
        data[1 * 16 + 0 * 4 + 0] = xtt; data[1 * 16 + 0 * 4 + 1] = xtx; data[1 * 16 + 0 * 4 + 2] = xty; data[1 * 16 + 0 * 4 + 3] = xtz;
        data[1 * 16 + 1 * 4 + 0] = xxt; data[1 * 16 + 1 * 4 + 1] = xxx; data[1 * 16 + 1 * 4 + 2] = xxy; data[1 * 16 + 1 * 4 + 3] = xxz;
        data[1 * 16 + 2 * 4 + 0] = xyt; data[1 * 16 + 2 * 4 + 1] = xyx; data[1 * 16 + 2 * 4 + 2] = xyy; data[1 * 16 + 2 * 4 + 3] = xyz;
        data[1 * 16 + 3 * 4 + 0] = xzt; data[1 * 16 + 3 * 4 + 1] = xzx; data[1 * 16 + 3 * 4 + 2] = xzy; data[1 * 16 + 3 * 4 + 3] = xzz;
        data[2 * 16 + 0 * 4 + 0] = ytt; data[2 * 16 + 0 * 4 + 1] = ytx; data[2 * 16 + 0 * 4 + 2] = yty; data[2 * 16 + 0 * 4 + 3] = ytz;
        data[2 * 16 + 1 * 4 + 0] = yxt; data[2 * 16 + 1 * 4 + 1] = yxx; data[2 * 16 + 1 * 4 + 2] = yxy; data[2 * 16 + 1 * 4 + 3] = yxz;
        data[2 * 16 + 2 * 4 + 0] = yyt; data[2 * 16 + 2 * 4 + 1] = yyx; data[2 * 16 + 2 * 4 + 2] = yyy; data[2 * 16 + 2 * 4 + 3] = yyz;
        data[2 * 16 + 3 * 4 + 0] = yzt; data[2 * 16 + 3 * 4 + 1] = yzx; data[2 * 16 + 3 * 4 + 2] = yzy; data[2 * 16 + 3 * 4 + 3] = yzz;
        data[3 * 16 + 0 * 4 + 0] = ztt; data[3 * 16 + 0 * 4 + 1] = ztx; data[3 * 16 + 0 * 4 + 2] = zty; data[3 * 16 + 0 * 4 + 3] = ztz;
        data[3 * 16 + 1 * 4 + 0] = zxt; data[3 * 16 + 1 * 4 + 1] = zxx; data[3 * 16 + 1 * 4 + 2] = zxy; data[3 * 16 + 1 * 4 + 3] = zxz;
        data[3 * 16 + 2 * 4 + 0] = zyt; data[3 * 16 + 2 * 4 + 1] = zyx; data[3 * 16 + 2 * 4 + 2] = zyy; data[3 * 16 + 2 * 4 + 3] = zyz;
        data[3 * 16 + 3 * 4 + 0] = zzt; data[3 * 16 + 3 * 4 + 1] = zzx; data[3 * 16 + 3 * 4 + 2] = zzy; data[3 * 16 + 3 * 4 + 3] = zzz;
    }
    
    // Accessors:
    double& operator[](const rank3Indices& index)
    { return data[index.i*16 + index.j*4 + index.k]; }
    const double& operator[](const rank3Indices& index) const
    { return data[index.i*16 + index.j*4 + index.k]; }
    
    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size,' ');
        std::cout << space << "   (" << Format(data[0 * 16 + 0 * 4 + 0], precision) << "," << Format(data[0 * 16 + 0 * 4 + 1], precision) << "," << Format(data[0 * 16 + 0 * 4 + 2], precision) << "," << Format(data[0 * 16 + 0 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[0 * 16 + 1 * 4 + 0], precision) << "," << Format(data[0 * 16 + 1 * 4 + 1], precision) << "," << Format(data[0 * 16 + 1 * 4 + 2], precision) << "," << Format(data[0 * 16 + 1 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[0 * 16 + 2 * 4 + 0], precision) << "," << Format(data[0 * 16 + 2 * 4 + 1], precision) << "," << Format(data[0 * 16 + 2 * 4 + 2], precision) << "," << Format(data[0 * 16 + 2 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[0 * 16 + 3 * 4 + 0], precision) << "," << Format(data[0 * 16 + 3 * 4 + 1], precision) << "," << Format(data[0 * 16 + 3 * 4 + 2], precision) << "," << Format(data[0 * 16 + 3 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (--------------------------------)\n";
        std::cout << space << "   (" << Format(data[1 * 16 + 0 * 4 + 0], precision) << "," << Format(data[1 * 16 + 0 * 4 + 1], precision) << "," << Format(data[1 * 16 + 0 * 4 + 2], precision) << "," << Format(data[1 * 16 + 0 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[1 * 16 + 1 * 4 + 0], precision) << "," << Format(data[1 * 16 + 1 * 4 + 1], precision) << "," << Format(data[1 * 16 + 1 * 4 + 2], precision) << "," << Format(data[1 * 16 + 1 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[1 * 16 + 2 * 4 + 0], precision) << "," << Format(data[1 * 16 + 2 * 4 + 1], precision) << "," << Format(data[1 * 16 + 2 * 4 + 2], precision) << "," << Format(data[1 * 16 + 2 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[1 * 16 + 3 * 4 + 0], precision) << "," << Format(data[1 * 16 + 3 * 4 + 1], precision) << "," << Format(data[1 * 16 + 3 * 4 + 2], precision) << "," << Format(data[1 * 16 + 3 * 4 + 3], precision) << ")" << "\n";
        std::cout << name  << " = (--------------------------------)\n";
        std::cout << space << "   (" << Format(data[2 * 16 + 0 * 4 + 0], precision) << "," << Format(data[2 * 16 + 0 * 4 + 1], precision) << "," << Format(data[2 * 16 + 0 * 4 + 2], precision) << "," << Format(data[2 * 16 + 0 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 16 + 1 * 4 + 0], precision) << "," << Format(data[2 * 16 + 1 * 4 + 1], precision) << "," << Format(data[2 * 16 + 1 * 4 + 2], precision) << "," << Format(data[2 * 16 + 1 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 16 + 2 * 4 + 0], precision) << "," << Format(data[2 * 16 + 2 * 4 + 1], precision) << "," << Format(data[2 * 16 + 2 * 4 + 2], precision) << "," << Format(data[2 * 16 + 2 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[2 * 16 + 3 * 4 + 0], precision) << "," << Format(data[2 * 16 + 3 * 4 + 1], precision) << "," << Format(data[2 * 16 + 3 * 4 + 2], precision) << "," << Format(data[2 * 16 + 3 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (--------------------------------)\n";
        std::cout << space << "   (" << Format(data[3 * 16 + 0 * 4 + 0], precision) << "," << Format(data[3 * 16 + 0 * 4 + 1], precision) << "," << Format(data[3 * 16 + 0 * 4 + 2], precision) << "," << Format(data[3 * 16 + 0 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[3 * 16 + 1 * 4 + 0], precision) << "," << Format(data[3 * 16 + 1 * 4 + 1], precision) << "," << Format(data[3 * 16 + 1 * 4 + 2], precision) << "," << Format(data[3 * 16 + 1 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[3 * 16 + 2 * 4 + 0], precision) << "," << Format(data[3 * 16 + 2 * 4 + 1], precision) << "," << Format(data[3 * 16 + 2 * 4 + 2], precision) << "," << Format(data[3 * 16 + 2 * 4 + 3], precision) << ")" << "\n";
        std::cout << space << "   (" << Format(data[3 * 16 + 3 * 4 + 0], precision) << "," << Format(data[3 * 16 + 3 * 4 + 1], precision) << "," << Format(data[3 * 16 + 3 * 4 + 2], precision) << "," << Format(data[3 * 16 + 3 * 4 + 3], precision) << ")" << "\n";
        if(newline) std::cout << "\n";
    }
    friend std::ostream& operator<<(std::ostream& os, const Tensor4x4x4& v);
};
inline std::ostream& operator<<(std::ostream& os, const Tensor4x4x4& v) 
{
    os << "(" << v[{0,0,0}] << "," << v[{0,0,1}] << "," << v[{0,0,2}] << "," << v[{0,0,3}]
       << "|" << v[{0,1,0}] << "," << v[{0,1,1}] << "," << v[{0,1,2}] << "," << v[{0,1,3}]
       << "|" << v[{0,2,0}] << "," << v[{0,2,1}] << "," << v[{0,2,2}] << "," << v[{0,2,3}]
       << "|" << v[{0,3,0}] << "," << v[{0,3,1}] << "," << v[{0,3,2}] << "," << v[{0,3,3}] << ")\n";
    os << "(" << v[{1,0,0}] << "," << v[{1,0,1}] << "," << v[{1,0,2}] << "," << v[{1,0,3}]
       << "|" << v[{1,1,0}] << "," << v[{1,1,1}] << "," << v[{1,1,2}] << "," << v[{1,1,3}]
       << "|" << v[{1,2,0}] << "," << v[{1,2,1}] << "," << v[{1,2,2}] << "," << v[{1,2,3}]
       << "|" << v[{1,3,0}] << "," << v[{1,3,1}] << "," << v[{1,3,2}] << "," << v[{1,3,3}] << ")\n";
    os << "(" << v[{2,0,0}] << "," << v[{2,0,1}] << "," << v[{2,0,2}] << "," << v[{2,0,3}]
       << "|" << v[{2,1,0}] << "," << v[{2,1,1}] << "," << v[{2,1,2}] << "," << v[{2,1,3}]
       << "|" << v[{2,2,0}] << "," << v[{2,2,1}] << "," << v[{2,2,2}] << "," << v[{2,2,3}]
       << "|" << v[{2,3,0}] << "," << v[{2,3,1}] << "," << v[{2,3,2}] << "," << v[{2,3,3}] << ")\n";
    os << "(" << v[{3,0,0}] << "," << v[{3,0,1}] << "," << v[{3,0,2}] << "," << v[{3,0,3}]
       << "|" << v[{3,1,0}] << "," << v[{3,1,1}] << "," << v[{3,1,2}] << "," << v[{3,1,3}]
       << "|" << v[{3,2,0}] << "," << v[{3,2,1}] << "," << v[{3,2,2}] << "," << v[{3,2,3}]
       << "|" << v[{3,3,0}] << "," << v[{3,3,1}] << "," << v[{3,3,2}] << "," << v[{3,3,3}] << ")";
    return os;
}





template<typename T>
struct UnstructuredMatrix
{
private:
    std::vector<T,AlignedArrayAllocator<T>> entries;
    std::vector<size_t,AlignedArrayAllocator<size_t>> indexes;

public:
    UnstructuredMatrix()
    { indexes.push_back(0); }

    void AddRow(T* data, size_t count)
    {
        for(size_t i=0; i<count; i++)
            entries.push_back(data[i]);
        indexes.push_back(indexes[indexes.size() - 1] + count);
    }
    void AddRow(std::vector<T> data)
    {
        for(size_t i=0; i<data.size(); i++)
            entries.push_back(data[i]);
        indexes.push_back(indexes[indexes.size() - 1] + data.size());
    }

    size_t Start(size_t row)
    { return indexes[row]; }
    size_t End(size_t row)
    { return indexes[row+1]; }

    T& operator[](const size_t index)
    { return entries[index]; }
    const T& operator[](const size_t index) const
    { return entries[index]; }
    T& operator[](const rank2Indices& index)
    { return entries[indexes[index.i] + index.j]; }
    const T& operator[](const rank2Indices& index) const
    { return entries[indexes[index.i] + index.j]; }

    void Print()
    {
        std::cout << "index: [entries]\n";
        for(size_t i=0; i<indexes.size()-1; i++)   // last entry outside for loop
        {
            std::cout << i << ": [";
            for(size_t j=indexes[i]; j<indexes[i+1] - 1; j++)  // last entry outside for loop
                std::cout << entries[j] << ", ";
            std::cout << entries[indexes[i+1] - 1];
            std::cout << "]\n";
        }
    }

    void PrintFlat()
    {
        std::cout << "entries: ";
        for(size_t i=0; i<entries.size(); i++)
            std::cout << entries[i] << ", ";
        std::cout << std::endl;
        std::cout << "indexes: ";
        for(size_t i=0; i<indexes.size(); i++)
            std::cout << indexes[i] << ", ";
        std::cout << std::endl;
    }
};
#endif //__INCLUDE_GUARD_TesorTypes_hh__