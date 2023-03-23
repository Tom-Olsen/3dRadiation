#include "Stencil.h"



// ------------------------------- Stencil -------------------------------
void Stencil::SetCoefficientCount()
{ nCoefficients = ((nOrder + 1) / 2) * ((nOrder + 1) / 2); }
void Stencil::AllocateBuffers()
{
    w.resize(nDir);
    cx.resize(nDir);
    cy.resize(nDir);
    cz.resize(nDir);
    theta.resize(nDir);
    phi.resize(nDir);
}
double Stencil::W(size_t d) const
{ return w[d]; }
double Stencil::Theta(size_t d) const
{ return theta[d]; }
double Stencil::Phi(size_t d) const
{ return phi[d]; }
double Stencil::Cx(size_t d) const
{ return cx[d]; }
double Stencil::Cy(size_t d) const
{ return cy[d]; }
double Stencil::Cz(size_t d) const
{ return cz[d]; }
Tensor3 Stencil::C(size_t d) const
{ return Tensor3(Cx(d), Cy(d), Cz(d)); }
// -----------------------------------------------------------------------



// ------------------------------ MyStencil ------------------------------
MyStencil::MyStencil(size_t nOrder)
{
    nTh = nOrder;
    nPh = 2 * nOrder;
    this->nOrder = nOrder;
    this->nDir = nTh * nPh;
    SetCoefficientCount();

    switch(nOrder)
    {
        case  3: dTheta = 0.88607712379261370; break;
        case  5: dTheta = 0.56708068878468730; break;
        case  7: dTheta = 0.41679707883494480; break;
        case  9: dTheta = 0.32944347754574150; break;
        case 11: dTheta = 0.27234940787623110; break;
        case 13: dTheta = 0.23211697811143664; break;
        case 15: dTheta = 0.20223903140287983; break;
        case 17: dTheta = 0.17917455393695525; break;
        case 19: dTheta = 0.16083171772074056; break;
        default:
            exit_on_error("MyStencil, invalid nOrder.");
    }

    // Initialize weights:
    w0 = 0;
    for(size_t d1=0; d1<nPh; d1++)
    for(size_t d0=0; d0<nTh; d0++)
        w0 += MySin(Theta(d0,d1));
    w0 = 4.0 * M_PI / w0;
}
// Overrides:
double MyStencil::W(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return W(d0,d1);
}
double MyStencil::Theta(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Theta(d0,d1);
}
double MyStencil::Phi(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Phi(d0,d1);
}
double MyStencil::Cx(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Cx(d0,d1);
}
double MyStencil::Cy(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Cy(d0,d1);
}
double MyStencil::Cz(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Cz(d0,d1);
}
Tensor3 MyStencil::C(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Tensor3(Cx(d0,d1), Cy(d0,d1), Cz(d0,d1));
}

// Indexing:
size_t MyStencil::Index(size_t d0, size_t d1) const
{ return d0 + d1 * nTh; }
double MyStencil::d0(double theta) const
{ return nTh / 2 + (theta - M_PI_2) / dTheta; }
double MyStencil::d1(double phi) const
{ return phi * nPh / (2.0 * M_PI) - 0.5; }

// Acces with two indices:
double MyStencil::W(size_t d0, size_t d1) const
{ return w0 * MySin(Theta(d0,d1)); }

double MyStencil::Theta(double d0, double d1) const
{ return M_PI_2 + (d0 - nTh / 2) * dTheta; }
double MyStencil::Phi(double d0, double d1) const
{ return 2.0 * M_PI * (d1 + 0.5) / nPh; }

double MyStencil::Cx(double d0, double d1) const
{ return MySin(Theta(d0,d1)) * MyCos(Phi(d0,d1)); }
double MyStencil::Cy(double d0, double d1) const
{ return MySin(Theta(d0,d1)) * MySin(Phi(d0,d1)); }
double MyStencil::Cz(double d0, double d1) const
{ return MyCos(Theta(d0,d1)); }
Tensor3 MyStencil::C(double d0, double d1) const
{ return Tensor3(Cx(d0,d1), Cy(d0,d1), Cz(d0,d1)); }
// ---------------------------------------------------------------------



// -------------------------- LebedevStencil ---------------------------
LebedevStencil::LebedevStencil(size_t nOrder)
{
    this->nOrder = nOrder;
    SetCoefficientCount();

    switch (nOrder)
    {
        case  3: nDir =   6; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil3" 
                 break;
        case  5: nDir =  14; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil5" 
                 break;
        case  7: nDir =  26; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil7" 
                 break;
        case  9: nDir =  38; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil9" 
                 break;
        case 11: nDir =  50; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil11"
                 break;
        case 13: nDir =  74; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil13"
                 break;
        case 15: nDir =  86; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil15"
                 break;
        case 17: nDir = 110; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil17"
                 break;
        case 19: nDir = 146; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil19"
                 break;
        case 21: nDir = 170; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil21"
                 break;
        case 23: nDir = 194; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil23"
                 break;
        case 25: nDir = 230; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil25"
                 break;
        case 27: nDir = 266; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil27"
                 break;
        case 29: nDir = 302; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil29"
                 break;
        case 31: nDir = 350; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil31"
                 break;
        default:
            exit_on_error("Invalid LebedevStencil nOrder. Must be odd and smaller 32.");
            break;
    }
}
// ---------------------------------------------------------------------



// ----------------------- GaussLegendreStencil ------------------------
GaussLegendreStencil::GaussLegendreStencil(size_t nOrder)
{
    this->nOrder = nOrder;
    this->nDir = nOrder * (2 * nOrder);
    SetCoefficientCount();

    switch (nOrder)
    {
        case  3: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil3" 
                 break;
        case  5: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil5" 
                 break;
        case  7: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil7" 
                 break;
        case  9: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil9" 
                 break;
        case 11: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil11"
                 break;
        case 13: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil13"
                 break;
        case 15: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil15"
                 break;
        case 17: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil17"
                 break;
        case 19: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil19"
                 break;
        case 21: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil21"
                 break;
        case 23: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil23"
                 break;
        case 25: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil25"
                 break;
        case 27: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil27"
                 break;
        case 29: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil29"
                 break;
        case 31: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil31"
                 break;
        default:
            exit_on_error("Invalid GaussLegendreStencil nOrder. Must be odd and smaller 32.");
            break;
    }
}
// ---------------------------------------------------------------------