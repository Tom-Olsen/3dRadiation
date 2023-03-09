#ifndef __INCLUDE_GUARD_Stencil_hh__
#define __INCLUDE_GUARD_Stencil_hh__
#include <glm/gtc/quaternion.hpp>   // Quaternions.
#include <fstream>                  // File input/output.
#include "Utility.hh"
#include "TensorTypes.hh"



// Spherical Coordinates convention:
// theta € [0,pi]
//   phi € [0,2pi]
// theta =  0 => North Pole
// theta = pi => South Pole
// (th,phi) = (pi/2,    0) => (x,y,z) = ( 1, 0, 0)
// (th,phi) = (pi/2, pi/2) => (x,y,z) = ( 0, 1, 0)
// (th,phi) = (pi/2,   pi) => (x,y,z) = (-1, 1, 0)
// (th,phi) = (pi/2,3pi/2) => (x,y,z) = ( 0,-1, 0)
// (th,phi) = (pi/2,  2pi) => (x,y,z) = ( 1, 0, 0)






// This class holds discretized velocities and the corresponding weights.
// The velocity distribution is mimicing that ob the Gauß-Legendre Quadrature
// as close as possible, while still being a uniform grid on the unit sphere.
// This makes the quadrature sligly inaccurate, which is a tradeoff for a faster velcoity interpolation.
// nCoefficients = 1+3+5+...+nTh is the number of spherical harmonics which can be integrated 'exactly'.
struct Stencil
{
    size_t nTh; // = order
    size_t nPh;
    size_t nThPh;
    size_t nCoefficients;
    double w0;
    double sigma = 1.0;
    double dTheta;

    Stencil(size_t nTh) : nTh(nTh), nPh(2 * nTh), nThPh(nTh * nPh)
    {
        if      (nTh ==  3) { dTheta = 0.88607712379261370; }
        else if (nTh ==  5) { dTheta = 0.56708068878468730; }
        else if (nTh ==  7) { dTheta = 0.41679707883494480; }
        else if (nTh ==  9) { dTheta = 0.32944347754574150; }
        else if (nTh == 11) { dTheta = 0.27234940787623110; }
        else if (nTh == 13) { dTheta = 0.23211697811143664; }
        else if (nTh == 15) { dTheta = 0.20223903140287983; }
        else if (nTh == 17) { dTheta = 0.17917455393695525; }
        else if (nTh == 19) { dTheta = 0.16083171772074056; }
        else if (nTh % 2 == 0) 
            exit_on_error("Even nTh not allowed.");
        else
            exit_on_error("Given stencil nTh not implemented (yet).");

        // Set number of spherical harmonics which can be integrated:
        nCoefficients = 0;
        for(size_t i=1; i<=nTh; i+=2)
            nCoefficients += i;

        // Initialize weights:
        w0 = 0;
        for(size_t d1=0; d1<nPh; d1++)
        for(size_t d0=0; d0<nTh; d0++)
            w0 += MySin<9>(Theta(d0,d1));
        w0 = 4.0 * M_PI / w0;
    }

    size_t Index(size_t d0, size_t d1) const
    { return d0 + d1 * nTh; }
    double d0(double theta) const
    { return nTh / 2 + (theta - M_PI_2) / dTheta; }
    double d1(double phi) const
    { return phi * nPh / (2.0 * M_PI) - 0.5; }

    double Theta(double d0, double d1) const
    { return M_PI_2 + (d0 - nTh / 2) * dTheta; }
    double Phi(double d0, double d1) const
    { return 2.0 * M_PI * (d1 + 0.5) / nPh; }

    double W(size_t d0, size_t d1) const
    { return w0 * MySin<9>(Theta(d0,d1)); }

    // Velocity vector.
    double Cx(double d0, double d1) const
    { return MySin<9>(Theta(d0,d1)) * MyCos<9>(Phi(d0,d1)); }
    double Cy(double d0, double d1) const
    { return MySin<9>(Theta(d0,d1)) * MySin<9>(Phi(d0,d1)); }
    double Cz(double d0, double d1) const
    { return MyCos<9>(Theta(d0,d1)); }

    Tensor3 Cxyz(double d0, double d1) const
    { return Tensor3(Cx(d0,d1), Cy(d0,d1), Cz(d0,d1)); }



    void Print()
    {
        std::cout << "  theta,  phi,    w,      x,      y,      z:\n";
        for(size_t d1=0; d1<nPh; d1++)
        for(size_t d0=0; d0<nTh; d0++)
            std::cout << "(" << Format(Theta(d0,d1)) << ", " <<Format(Phi(d0,d1)) << ", " << Format(W(d0,d1)) <<
                        ", " << Format(Cx(d0,d1)) << ", " << Format(Cy(d0,d1)) << ", " << Format(Cz(d0,d1)) <<")\n";
        std::cout << std::endl;
    }

    void WriteToCsv()
    {
        std::ofstream file("output/stencil.csv");
        file << "#x,y,z,w\n";
        for(size_t d1=0; d1<nPh; d1++)
        for(size_t d0=0; d0<nTh; d0++)
            file << Format(Cx(d0,d1)) << ", " << Format(Cy(d0,d1)) << ", " << Format(Cz(d0,d1)) << "," << Format(W(d0,d1)) << "\n";
        file.close();
    }
};



// This class holds discretized velocities and the corresponding weights.
// The distribution of the velocity vectors is given by the spherical Lebedev quadarture.
// Exact Spherical Harmonic integration:
// LebedevStencil3: Y00, Y1m1,Y11,Y1p1
// LebedevStencil5: Y00, Y1m1,Y11,Y1p1, Y2m2,Y2m1,Y20,Y2p1,Y2p2
// LebedevStencil7: Y00, Y1m1,Y11,Y1p1, Y2m2,Y2m1,Y20,Y2p1,Y2p2, Y3m3,Y3m2,Y3m1,Y30,Y3p1,Y3p2,Y3p3
// nCoefficients = 1+3+5+...+n is the number of spherical harmonics which can be integrated exactly,
// where n in the number of the choosen LebedevStencil.
struct LebedevStencil
{
public:
    size_t nDir;
    size_t nOrder;
    size_t nCoefficients;
protected:
    double* w;
    double* cx;
    double* cy;
    double* cz;

public:
    double Theta(double x, double y, double z)
    { return MyAtan2(sqrt(x * x + y * y), z); }
    double Theta(size_t d)
    { return Theta(Cx(d), Cy(d), Cz(d)); }
    double Phi(double x, double y, double z)
    {
        double phi = MyAtan2(y, x);
        return (phi < 0) ? phi + 2 * M_PI : phi;
    }
    double Phi(size_t d)
    { return Phi(Cx(d), Cy(d), Cz(d)); }

    double W(size_t d) const
    { return w[d]; }

    // Velocity vector.
    double Cx(size_t d) const
    { return cx[d]; }
    double Cy(size_t d) const
    { return cy[d]; }
    double Cz(size_t d) const
    { return cz[d]; }
    Tensor3 Cxyz(size_t d) const
    { return Tensor3(cx[d], cy[d], cz[d]); }
    


    void Print()
    {
        std::cout << "  theta,  phi,    w,      x,      y,      z:\n";
        for(size_t d=0; d<nDir; d++)
            std::cout << "(" << Format(Theta(d)) << ", " <<Format(Phi(d)) << ", " << Format(W(d)) <<
                        ", " << Format(Cx(d)) << ", " << Format(Cy(d)) << ", " << Format(Cz(d)) <<")\n";
        std::cout << std::endl;
    }

    void WriteToCsv()
    {
        std::ofstream file("output/lebedevStencil.csv");
        file << "#x,y,z,w\n";
        for(size_t d=0; d<nDir; d++)
            file << Format(Cx(d)) << ", " << Format(Cy(d)) << ", " << Format(Cz(d)) << "," << Format(W(d)) << "\n";
        file.close();
    }
};

struct LebedevStencil3 : LebedevStencil
{
    LebedevStencil3()
    {
        this->nDir = 6;
        this->nOrder = 3;
        this->nCoefficients = 1 + 3;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil3"
    }
};
struct LebedevStencil5 : LebedevStencil
{
    LebedevStencil5()
    {
        this->nDir = 14;
        this->nOrder = 5;
        this->nCoefficients = 1 + 3 + 5;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil5"
    }
};
struct LebedevStencil7 : LebedevStencil
{
    LebedevStencil7()
    {
        this->nDir = 26;
        this->nOrder = 7;
        this->nCoefficients = 1 + 3 + 5 + 7;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil7"
    }
};
struct LebedevStencil9 : LebedevStencil
{
    LebedevStencil9()
    {
        this->nDir = 38;
        this->nOrder = 9;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil9"
    }
};
struct LebedevStencil11 : LebedevStencil
{
    LebedevStencil11()
    {
        this->nDir = 50;
        this->nOrder = 11;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil11"
    }
};
struct LebedevStencil13 : LebedevStencil
{
    LebedevStencil13()
    {
        this->nDir = 74;
        this->nOrder = 13;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11 + 13;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil13"
    }
};
struct LebedevStencil15 : LebedevStencil
{
    LebedevStencil15()
    {
        this->nDir = 86;
        this->nOrder = 15;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11 + 13 + 15;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "LebedevStencil/LebedevStencil15"
    }
};



// This type of stencil is not used in the 'general relativistic Lattice Boltzmann Method for radiative transport' code.
// It is only here for reference, if someone might want to test something.
// The 'Stencil' class above is an approximation of this stencil class.
// The approximation allows for a uniform grid which makes interpolations in velocity space easier.
struct GaussLegendreStencil
{
public:
    size_t nDir;
    size_t nOrder;
    size_t nCoefficients;
protected:
    double* w;
    double* cx;
    double* cy;
    double* cz;

public:
    double Theta(double x, double y, double z)
    { return MyAtan2(sqrt(x * x + y * y), z); }
    double Theta(size_t d)
    { return Theta(Cx(d), Cy(d), Cz(d)); }
    double Phi(double x, double y, double z)
    {
        double phi = MyAtan2(y, x);
        return (phi < 0) ? phi + 2 * M_PI : phi;
    }
    double Phi(size_t d)
    { return Phi(Cx(d), Cy(d), Cz(d)); }

    double W(size_t d) const
    { return w[d]; }

    // Velocity vector.
    double Cx(size_t d) const
    { return cx[d]; }
    double Cy(size_t d) const
    { return cy[d]; }
    double Cz(size_t d) const
    { return cz[d]; }
    Tensor3 Cxyz(size_t d) const
    { return Tensor3(cx[d], cy[d], cz[d]); }
    


    void Print()
    {
        std::cout << "  theta,  phi,    w,      x,      y,      z:\n";
        for(size_t d=0; d<nDir; d++)
            std::cout << "(" << Format(Theta(d)) << ", " <<Format(Phi(d)) << ", " << Format(W(d)) <<
                        ", " << Format(Cx(d)) << ", " << Format(Cy(d)) << ", " << Format(Cz(d)) <<")\n";
        std::cout << std::endl;
    }

    void WriteToCsv()
    {
        std::ofstream file("output/gaussLegendreStencil.csv");
        file << "#x,y,z,w\n";
        for(size_t d=0; d<nDir; d++)
            file << Format(Cx(d)) << ", " << Format(Cy(d)) << ", " << Format(Cz(d)) << "," << Format(W(d)) << "\n";
        file.close();
    }
};
struct GaussLegendreStencil3 : GaussLegendreStencil
{
    GaussLegendreStencil3()
    {
        this->nDir = 18;
        this->nOrder = 3;
        this->nCoefficients = 1 + 3;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil3"
    }
};
struct GaussLegendreStencil5 : GaussLegendreStencil
{
    GaussLegendreStencil5()
    {
        this->nDir = 50;
        this->nOrder = 5;
        this->nCoefficients = 1 + 3 + 5;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil5"
    }
};
struct GaussLegendreStencil7 : GaussLegendreStencil
{
    GaussLegendreStencil7()
    {
        this->nDir = 98;
        this->nOrder = 7;
        this->nCoefficients = 1 + 3 + 5 + 7;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil7"
    }
};
struct GaussLegendreStencil9 : GaussLegendreStencil
{
    GaussLegendreStencil9()
    {
        this->nDir = 162;
        this->nOrder = 9;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil9"
    }
};
struct GaussLegendreStencil11 : GaussLegendreStencil
{
    GaussLegendreStencil11()
    {
        this->nDir = 242;
        this->nOrder = 11;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil11"
    }
};
struct GaussLegendreStencil13 : GaussLegendreStencil
{
    GaussLegendreStencil13()
    {
        this->nDir = 338;
        this->nOrder = 13;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11 + 13;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil13"
    }
};
struct GaussLegendreStencil15 : GaussLegendreStencil
{
    GaussLegendreStencil15()
    {
        this->nDir = 450;
        this->nOrder = 15;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11 + 13 + 15;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil15"
    }
};
struct GaussLegendreStencil17 : GaussLegendreStencil
{
    GaussLegendreStencil17()
    {
        this->nDir = 578;
        this->nOrder = 17;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11 + 13 + 15 + 17;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil17"
    }
};
struct GaussLegendreStencil19 : GaussLegendreStencil
{
    GaussLegendreStencil19()
    {
        this->nDir = 722;
        this->nOrder = 19;
        this->nCoefficients = 1 + 3 + 5 + 7 + 9 + 11 + 13 + 15 + 17 + 19;
        this->w = new double[nDir];
        this->cx = new double[nDir];
        this->cy = new double[nDir];
        this->cz = new double[nDir];
        #include "GaussLegendreStencil/GaussLegendreStencil19"
    }
};
#endif //__INCLUDE_GUARD_Stencil_hh__