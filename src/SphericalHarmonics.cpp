#include "SphericalHarmonics.h"

// ------------------------------ Spherical Harmonics ------------------------------
// l = 0:
double SphericalHarmonics::Y00(double x, double y, double z)
{ return 0.5 * sqrt(1.0 / M_PI); }

double SphericalHarmonics::Y00(const Tensor3& xyz)
{ return 0.5 * sqrt(1.0 / M_PI); }

// l = 1:
double SphericalHarmonics::Y1m1(double x, double y, double z)
{ return sqrt(3.0 / (4.0 * M_PI)) * y; }
double SphericalHarmonics::Y10 (double x, double y, double z)
{ return sqrt(3.0 / (4.0 * M_PI)) * z; }
double SphericalHarmonics::Y1p1(double x, double y, double z)
{ return sqrt(3.0 / (4.0 * M_PI)) * x; }

double SphericalHarmonics::Y1m1(const Tensor3& xyz)
{ return Y1m1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y10 (const Tensor3& xyz)
{ return Y10 (xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y1p1(const Tensor3& xyz)
{ return Y1p1(xyz[1], xyz[2], xyz[3]); }



// l = 2:
double SphericalHarmonics::Y2m2(double x, double y, double z)
{ return 0.5  * sqrt(15.0 / M_PI) * x * y; }
double SphericalHarmonics::Y2m1(double x, double y, double z)
{ return 0.5  * sqrt(15.0 / M_PI) * y * z; }
double SphericalHarmonics::Y20 (double x, double y, double z)
{ return 0.25 * sqrt( 5.0 / M_PI) * (2.0 * z * z - x * x - y * y); }
double SphericalHarmonics::Y2p1(double x, double y, double z)
{ return 0.5  * sqrt(15.0 / M_PI) * z * x; }
double SphericalHarmonics::Y2p2(double x, double y, double z)
{ return 0.25 * sqrt(15.0 / M_PI) * (x * x - y * y); }

double SphericalHarmonics::Y2m2(const Tensor3& xyz)
{ return Y2m2(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y2m1(const Tensor3& xyz)
{ return Y2m1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y20 (const Tensor3& xyz)
{ return Y20 (xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y2p1(const Tensor3& xyz)
{ return Y2p1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y2p2(const Tensor3& xyz)
{ return Y2p2(xyz[1], xyz[2], xyz[3]); }



// l = 3:
double SphericalHarmonics::Y3m3(double x, double y, double z)
{ return 0.25 * sqrt( 35.0 / (2.0 * M_PI)) * (3.0 * x * x - y * y) * y; }
double SphericalHarmonics::Y3m2(double x, double y, double z)
{ return 0.5  * sqrt(105.0 / M_PI) * x * y * z; }
double SphericalHarmonics::Y3m1(double x, double y, double z)
{ return 0.25 * sqrt( 21.0 / (2.0 * M_PI)) * (4.0 * z * z - x * x - y * y) * y; }
double SphericalHarmonics::Y30 (double x, double y, double z)
{ return 0.25 * sqrt(  7.0 / M_PI) * (2.0 * z * z - 3.0 * x * x - 3.0 * y * y) * z; }
double SphericalHarmonics::Y3p1(double x, double y, double z)
{ return 0.25 * sqrt( 21.0 / (2.0 * M_PI)) * (4.0 * z * z - x * x - y * y) * x; }
double SphericalHarmonics::Y3p2(double x, double y, double z)
{ return 0.25 * sqrt(105.0 / M_PI) * (x * x - y * y) * z; }
double SphericalHarmonics::Y3p3(double x, double y, double z)
{ return 0.25 * sqrt( 35.0/ (2.0 * M_PI)) * (x * x - 3.0 * y * y) * x; }

double SphericalHarmonics::Y3m3(const Tensor3& xyz)
{ return Y3m3(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y3m2(const Tensor3& xyz)
{ return Y3m2(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y3m1(const Tensor3& xyz)
{ return Y3m1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y30 (const Tensor3& xyz)
{ return Y30 (xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y3p1(const Tensor3& xyz)
{ return Y3p1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y3p2(const Tensor3& xyz)
{ return Y3p2(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y3p3(const Tensor3& xyz)
{ return Y3p3(xyz[1], xyz[2], xyz[3]); }



// l = 4:
double SphericalHarmonics::Y4m4(double x, double y, double z)
{ return 0.75 * sqrt(35.0 / M_PI) * (x * x - y * y) * x * y; }
double SphericalHarmonics::Y4m3(double x, double y, double z)
{ return 0.75 * sqrt(35.0 / (2.0 * M_PI)) * (3.0 * x * x - y * y) * y * z; }
double SphericalHarmonics::Y4m2(double x, double y, double z)
{ return 0.75 * sqrt(5.0 / M_PI) * (7.0 * z * z - 1.0) * x * y; }
double SphericalHarmonics::Y4m1(double x, double y, double z)
{ return 0.75 * sqrt(5.0 / (2.0 * M_PI)) * (7.0 * z * z - 3.0) * y * z; }
double SphericalHarmonics::Y40 (double x, double y, double z)
{ return 0.1875 * sqrt(1.0 / M_PI) * (35.0 * z * z * z * z - 30.0 * z * z + 3.0); }
double SphericalHarmonics::Y4p1(double x, double y, double z)
{ return 0.75 * sqrt(5.0 / (2.0 * M_PI)) * (7.0 * z * z - 3.0) * x * z; }
double SphericalHarmonics::Y4p2(double x, double y, double z)
{ return 0.375 * sqrt(5.0 / M_PI) * (x * x - y * y) * (7.0 * z * z - 1.0); }
double SphericalHarmonics::Y4p3(double x, double y, double z)
{ return 0.75 * sqrt(35.0 / (2.0 * M_PI)) * (x * x - 3.0 * y * y) * x * z; }
double SphericalHarmonics::Y4p4(double x, double y, double z)
{ return 0.1875 * sqrt(35.0 / M_PI) * ((x * x - 3.0 * y * y) * x * x - (3.0 * x * x - y * y) * y * y); }

double SphericalHarmonics::Y4m4(const Tensor3& xyz)
{ return Y4m4(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4m3(const Tensor3& xyz)
{ return Y4m3(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4m2(const Tensor3& xyz)
{ return Y4m2(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4m1(const Tensor3& xyz)
{ return Y4m1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y40 (const Tensor3& xyz)
{ return Y40 (xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4p1(const Tensor3& xyz)
{ return Y4p1(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4p2(const Tensor3& xyz)
{ return Y4p2(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4p3(const Tensor3& xyz)
{ return Y4p3(xyz[1], xyz[2], xyz[3]); }
double SphericalHarmonics::Y4p4(const Tensor3& xyz)
{ return Y4p4(xyz[1], xyz[2], xyz[3]); }



std::vector<double (*)(double, double, double)> SphericalHarmonics::Harmonics =
{ &SphericalHarmonics::Y00, 
  &SphericalHarmonics::Y1m1, &SphericalHarmonics::Y10 , &SphericalHarmonics::Y1p1,
  &SphericalHarmonics::Y2m2, &SphericalHarmonics::Y2m1, &SphericalHarmonics::Y20 , &SphericalHarmonics::Y2p1, &SphericalHarmonics::Y2p2,
  &SphericalHarmonics::Y3m3, &SphericalHarmonics::Y3m2, &SphericalHarmonics::Y3m1, &SphericalHarmonics::Y30 , &SphericalHarmonics::Y3p1, &SphericalHarmonics::Y3p2, &SphericalHarmonics::Y3p3,
  &SphericalHarmonics::Y4m4, &SphericalHarmonics::Y4m3, &SphericalHarmonics::Y4m2, &SphericalHarmonics::Y4m1, &SphericalHarmonics::Y40 , &SphericalHarmonics::Y4p1, &SphericalHarmonics::Y4p2, &SphericalHarmonics::Y4p3, &SphericalHarmonics::Y4p4};
// ---------------------------------------------------------------------------------



// ------------------------------ Spherical Harmonics Expansion 3 ------------------------------
std::vector<double> SphericalHarmonics::GetCoefficients(const Stencil& stencil, const double* data, size_t nCoefficients)
{
    std::vector<double> coefficients(nCoefficients);
    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        size_t d = stencil.Index(d0,d1);
        double x = stencil.Cx(d0,d1);
        double y = stencil.Cy(d0,d1);
        double z = stencil.Cz(d0,d1);
        double c = data[d] * stencil.W(d0,d1);

        for(size_t i=0; i<nCoefficients; i++)
            coefficients[i] += c * Harmonics[i](x,y,z);
    }
    return coefficients;
}
std::vector<double> SphericalHarmonics::GetCoefficients(const LebedevStencil& stencil, const double* data, size_t nCoefficients)
{
    std::vector<double> coefficients(nCoefficients);
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double x = stencil.Cx(d);
        double y = stencil.Cy(d);
        double z = stencil.Cz(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<nCoefficients; i++)
            coefficients[i] += c * Harmonics[i](x,y,z);
    }
    return coefficients;
}
std::vector<double> SphericalHarmonics::GetCoefficients(const GaussLegendreStencil& stencil, const double* data, size_t nCoefficients)
{
    std::vector<double> coefficients(nCoefficients);
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double x = stencil.Cx(d);
        double y = stencil.Cy(d);
        double z = stencil.Cz(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<nCoefficients; i++)
            coefficients[i] += c * Harmonics[i](x,y,z);
    }
    return coefficients;
}

void SphericalHarmonics::GetCoefficients(const Stencil& stencil, const double* data, size_t nCoefficients, double* coefficients)
{
    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        size_t d = stencil.Index(d0,d1);
        double x = stencil.Cx(d0,d1);
        double y = stencil.Cy(d0,d1);
        double z = stencil.Cz(d0,d1);
        double c = data[d] * stencil.W(d0,d1);

        for(size_t i=0; i<nCoefficients; i++)
            coefficients[i] += c * Harmonics[i](x,y,z);
    }
}
void SphericalHarmonics::GetCoefficients(const LebedevStencil& stencil, const double* data, size_t nCoefficients, double* coefficients)
{
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double x = stencil.Cx(d);
        double y = stencil.Cy(d);
        double z = stencil.Cz(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<nCoefficients; i++)
            coefficients[i] += c * Harmonics[i](x,y,z);
    }
}
void SphericalHarmonics::GetCoefficients(const GaussLegendreStencil& stencil, const double* data, size_t nCoefficients, double* coefficients)
{
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double x = stencil.Cx(d);
        double y = stencil.Cy(d);
        double z = stencil.Cz(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<nCoefficients; i++)
            coefficients[i] += c * Harmonics[i](x,y,z);
    }
}



double SphericalHarmonics::GetValue(Tensor3 direction, const std::vector<double>& coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Harmonics[i](direction[1], direction[2], direction[3]);

    return result;
}
double SphericalHarmonics::GetValue(Tensor3 direction, double* coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Harmonics[i](direction[1], direction[2], direction[3]);

    return result;
}
double SphericalHarmonics::GetValue(double theta, double phi, const std::vector<double>& coefficients, size_t nCoefficients)
{
    double x = MySin(theta) * MyCos(phi);
    double y = MySin(theta) * MySin(phi);
    double z = MyCos(theta);
    
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Harmonics[i](x,y,z);

    return result;
}
double SphericalHarmonics::GetValue(double theta, double phi, double* coefficients, size_t nCoefficients)
{
    double x = MySin(theta) * MyCos(phi);
    double y = MySin(theta) * MySin(phi);
    double z = MyCos(theta);
    
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Harmonics[i](x,y,z);

    return result;
}
// ---------------------------------------------------------------------------------------------