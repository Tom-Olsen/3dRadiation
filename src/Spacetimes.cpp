#include "Spacetimes.h"

// ------------------------------ Minkowski ------------------------------
Minkowski::Minkowski(Grid& grid_, double m_, double a_) : Metric(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

bool Minkowski::InsideBH(const Coord& xyz)
{
    return false;
}
std::string Minkowski::Name()
{
    return "Minkowski";
}

Tensor4x4 Minkowski::MetricFunction(const Coord& xyz)
{
    return Tensor4x4(-1,0,0,0 ,0,1,0,0, 0,0,1,0, 0,0,0,1);
}
// -----------------------------------------------------------------------



// ------------------------------ SchwarzSchild ------------------------------
SchwarzSchild::SchwarzSchild(Grid& grid_, double m_, double a_) : Metric(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

bool SchwarzSchild::InsideBH(const Coord& xyz)
{
    // Buffer zone must be bigger here or geodesic equations dont converge.
    return xyz.EuklNormSquared() <= 2.5 * 2.5 * this->m * this->m;
}
std::string SchwarzSchild::Name()
{
    return "SchwarzSchild";
}

Tensor4x4 SchwarzSchild::MetricFunction(const Coord& xyz)
{
    double rs = 2.0 * this->m;
    double r2 = xyz.EuklNormSquared();
    double r = sqrt(r2);
    if(r > 2.0*this->m)
    {
        Tensor4x4 g_ll(0.0);
        g_ll[{0,0}] = -1.0 + rs / r;
        g_ll[{1,1}] = xyz[1]*xyz[1]/(r2-rs*r) + xyz[2]*xyz[2]/r2        + xyz[3]*xyz[3]/r2;
        g_ll[{2,2}] = xyz[1]*xyz[1]/r2        + xyz[2]*xyz[2]/(r2-rs*r) + xyz[3]*xyz[3]/r2;
        g_ll[{3,3}] = xyz[1]*xyz[1]/r2        + xyz[2]*xyz[2]/r2        + xyz[3]*xyz[3]/(r2-rs*r);
        g_ll[{2,1}] = g_ll[{1,2}] = rs*xyz[1]*xyz[2] / (r2*(r-rs));
        g_ll[{3,1}] = g_ll[{1,3}] = rs*xyz[1]*xyz[3] / (r2*(r-rs));
        g_ll[{3,2}] = g_ll[{2,3}] = rs*xyz[2]*xyz[3] / (r2*(r-rs));
        return g_ll;
    }
    else
        return Tensor4x4(-1,0,0,0 ,0,1,0,0, 0,0,1,0, 0,0,0,1);
}
// ---------------------------------------------------------------------------



// ------------------------------ KerrSchild ------------------------------
KerrSchild::KerrSchild(Grid& grid_, double m_, double a_) : Metric(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

bool KerrSchild::InsideBH(const Coord& xyz)
{
    return xyz.EuklNormSquared() <= 2.5 * 2.5 * this->m * this->m;
}
std::string KerrSchild::Name()
{
    return "KerrSchild";
}

Tensor4x4 KerrSchild::MetricFunction(const Coord& xyz)
{
    double a2 = this->a * this->a;
    double R2 = xyz.EuklNormSquared();
    if(R2 > a2)
    {
        double r2 = R2 - a2;
        double r  = sqrt(r2);
        double rho2 = r2 + a2;
        double H = this->m * IntegerPow<3>(r) / (IntegerPow<4>(r) + a2 * xyz[3] * xyz[3]);
        Tensor4 l_l( 1, (r * xyz[1] + this->a * xyz[2]) / rho2, (r * xyz[2] - this->a * xyz[1]) / rho2, xyz[3]/r);
        Tensor4 l_u(-1, (r * xyz[1] + this->a * xyz[2]) / rho2, (r * xyz[2] - this->a * xyz[1]) / rho2, xyz[3]/r);
        Tensor4x4 g_ll;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                if (i == 0 && j == 0)
                    g_ll[{i,j}] = -1.0 + 2.0 * H * l_l[i] * l_l[j];
                else
                    g_ll[{i,j}] = (i==j) + 2.0 * H * l_l[i] * l_l[j];
            }
        return g_ll;
    }
    else
        return Tensor4x4(-1,0,0,0 ,0,1,0,0, 0,0,1,0, 0,0,0,1);
}
// ------------------------------------------------------------------------