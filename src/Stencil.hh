#ifndef __INCLUDE_GUARD_Stencil_hh__
#define __INCLUDE_GUARD_Stencil_hh__
#include <glm/gtc/quaternion.hpp>   // Quaternions.
#include <glm/gtx/quaternion.hpp>   // Quaternions.
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
// The distribution of the velocity vectors is given by a spherical coordinate grid.
struct Stencil
{
    int nTh;
    int nPh;
    int nThPh;
    double w0;
    double sigma = 1.0; // more research is needed

    Stencil(int nTh, int nPh) : nTh(nTh), nPh(nPh)
    {// Normalize weights:
        nThPh = nTh * nPh;
        w0 = 0;
        for(int d1=0; d1<nPh; d1++)
        for(int d0=0; d0<nTh; d0++)
            w0 += sin(Theta(d0,d1));
        w0 = 4.0 * M_PI / w0;
    }

    int Index(int d0, int d1) const
    { return d0 + d1 * nTh; }
    double d0(double theta) const
    { return theta * nTh / M_PI - 0.5; }
    double d1(double phi) const
    { return phi * nPh / (2.0 * M_PI) - 0.5; }

    double Theta(double d0, double d1) const
    { return M_PI * (d0 + 0.5) / nTh; }
    double Theta(double x, double y, double z)
    { return atan2(sqrt(x * x + y * y), z); }

    double Phi(double d0, double d1) const
    { return 2.0 * M_PI * (d1 + 0.5) / nPh; }
    double Phi(double x, double y, double z)
    {
        double phi = atan2(y, x);
        return (phi < 0) ? phi + 2 * M_PI : phi;
    }

    double W(int d0, int d1) const
    { return w0 * sin(Theta(d0,d1)); }

    // Velocity vector.
    double Cx(double d0, double d1) const
    { return sin(Theta(d0,d1)) * cos(Phi(d0,d1)); }
    double Cy(double d0, double d1) const
    { return sin(Theta(d0,d1)) * sin(Phi(d0,d1)); }
    double Cz(double d0, double d1) const
    { return cos(Theta(d0,d1)); }

    Tensor3 Cxyz(double d0, double d1) const
    { return Tensor3(Cx(d0,d1), Cy(d0,d1), Cz(d0,d1)); }



    void Print()
    {
        std::cout << "  theta,  phi,    w,      x,      y,      z:\n";
        for(int d1=0; d1<nPh; d1++)
        for(int d0=0; d0<nTh; d0++)
            std::cout << "(" << Format(Theta(d0,d1),3) << ", " <<Format(Phi(d0,d1),3) << ", " << Format(W(d0,d1),3) <<
                        ", " << Format(Cx(d0,d1),3) << ", " << Format(Cy(d0,d1),3) << ", " << Format(Cz(d0,d1),3) <<")\n";
        std::cout << std::endl;
    }
};



// This class holds discretized velocities and the corresponding weights.
// The distribution of the velocity vectors is given by the spherical Lebedev quadarture.
// Accuracy of Spherical Harmonic integration:
// LebedevStencil3: Y00, Y1m1,Y11,Y1p1
// LebedevStencil5: Y00, Y1m1,Y11,Y1p1, Y2m2,Y2m1,Y20,Y2p1,Y2p2
// LebedevStencil7: Y00, Y1m1,Y11,Y1p1, Y2m2,Y2m1,Y20,Y2p1,Y2p2, Y3m3,Y3m2,Y3m1,Y30,Y3p1,Y3p2,Y3p3
struct LebedevStencil
{
public:
    int nDir;
    int nOrder;
    int nCoefficients;
protected:
    double* w;
    double* cx;
    double* cy;
    double* cz;

public:
    double Theta(double x, double y, double z)
    { return atan2(sqrt(x * x + y * y), z); }
    double Theta(int d)
    { return Theta(Cx(d), Cy(d), Cz(d)); }
    double Phi(double x, double y, double z)
    {
        double phi = atan2(y, x);
        return (phi < 0) ? phi + 2 * M_PI : phi;
    }
    double Phi(int d)
    { return Phi(Cx(d), Cy(d), Cz(d)); }

    double W(int d) const
    { return w[d]; }

    // Velocity vector.
    double Cx(int d) const
    { return cx[d]; }
    double Cy(int d) const
    { return cy[d]; }
    double Cz(int d) const
    { return cz[d]; }
    Tensor3 Cxyz(int d) const
    { return Tensor3(cx[d], cy[d], cz[d]); }
    


    void Print()
    {
        std::cout << "  theta,  phi,    w,      x,      y,      z:\n";
        for(int d=0; d<nDir; d++)
            std::cout << "(" << Format(Theta(d),3) << ", " <<Format(Phi(d),3) << ", " << Format(W(d),3) <<
                        ", " << Format(Cx(d),3) << ", " << Format(Cy(d),3) << ", " << Format(Cz(d),3) <<")\n";
        std::cout << std::endl;
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
#endif //__INCLUDE_GUARD_Stencil_hh__