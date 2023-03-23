#include <iostream>
#include "Utility.hh"
#include "TensorTypes.hh"
#include "glm/glm/gtc/quaternion.hpp"



struct Complex
{
    double real;
    double imag;
    
    Complex() : real(0), imag(0) {}
    Complex(double r, double i) : real(r), imag(i) {}
    static Complex Polar(double length, double angle)
    { return Complex(length * MyCos(angle), length * MySin(angle)); }

    Complex operator+(const Complex& other) const
    { return Complex(real + other.real, imag + other.imag); }
    
    Complex operator-(const Complex& other) const
    { return Complex(real - other.real, imag - other.imag); }
    
    Complex operator*(const Complex& other) const
    {
        double r = real * other.real - imag * other.imag;
        double i = real * other.imag + imag * other.real;
        return Complex(r, i);
    }
    Complex& operator +=(const Complex& other)
    {
        real += other.real;
        imag += other.imag;
        return *this;
    }
    
    Complex operator/(const Complex& other) const
    {
        double denom = other.real * other.real + other.imag * other.imag;
        double r = (real * other.real + imag * other.imag) / denom;
        double i = (imag * other.real - real * other.imag) / denom;
        return Complex(r, i);
    }
    
    bool operator==(const Complex& other) const
    { return real == other.real && imag == other.imag; }
    
    bool operator!=(const Complex& other) const
    { return !(*this == other); }
    
    double abs() const
    { return sqrt(real * real + imag * imag); }

    double angle() const
    { return MyAtan2(imag, real); }
    
    Complex conj() const
    { return Complex(real, -imag); }

    friend std::ostream& operator <<(std::ostream& os, const Complex& c);
};
inline std::ostream& operator <<(std::ostream& os, const Complex& c)
{
    os << c.real << " + " << c.imag << "i";
    return os;
}



struct Complex2
{
    Complex x;
    Complex y;

    Complex2() : x(Complex()), y(Complex()) {}
    Complex2(Complex x, Complex y) : x(x), y(y) {}

    Complex& operator[](const int index)
    {
        switch (index)
        {
            case 0: return x;
            case 1: return y;
        }
        exit_on_error("Invalid index in Complex2[index].");
        return x;
    }
    const Complex& operator[](const int index) const
    {
        switch (index)
        {
            case 0: return x;
            case 1: return y;
        }
        exit_on_error("Invalid index in Complex2[index].");
        return x;
    }
};



struct Complex2x2
{
    Complex xx;
    Complex xy;
    Complex yx;
    Complex yy;

    Complex2x2(Complex xx, Complex xy, Complex yx, Complex yy) : xx(xx), xy(xy), yx(yx), yy(yy) {}

    Complex& operator[](const rank2Indices& index)
    {
        switch (index.i*2 + index.j)
        {
            case 0: return xx;
            case 1: return xy;
            case 2: return yx;
            case 3: return yy;
        }
        exit_on_error("Invalid index in Complex2x2[index].");
        return xx;
    }
    const Complex& operator[](const rank2Indices& index) const
    {
        switch (index.i + 2*index.j)
        {
            case 0: return xx;
            case 1: return xy;
            case 2: return yx;
            case 3: return yy;
        }
        exit_on_error("Invalid index in Complex2x2[index].");
        return xx;
    }
};



inline Complex2 operator*(Complex2x2 mat, Complex2 vec)
{
    Complex2 result;
    for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
        result[i] += mat[{i,j}] * vec[j];
    return result;
}



inline void Rotate(double& theta, double& phi, const glm::quat& q)
{// https://stla.github.io/stlapblog/posts/RotationSphericalCoordinates.html
    Complex2 psi(Complex(MyCos(theta/2.0),0.0), Complex::Polar(MySin(theta/2.0),phi));
    Complex2x2 R
    ( Complex( q.w,-q.z), Complex(-q.y,-q.x),
      Complex( q.y,-q.x), Complex( q.w, q.z) );

    psi = R * psi;
    
    theta = 2.0 * MyAtan2(psi.y.abs(), psi.x.abs());
    phi = fmod(psi.y.angle() - psi.x.angle() + 2.0 * M_PI, 2.0 * M_PI);
}