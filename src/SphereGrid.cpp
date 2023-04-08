#include "SphereGrid.h"



SphereGrid::SphereGrid(size_t nTh, size_t nPh) : nTh(nTh), nPh(nPh), nDir(nTh*nPh) {}

double SphereGrid::Theta(size_t i) const
// { return M_PI * (i + 0.5) / (double)nTh; }
{ return M_PI * i / (nTh - 1.0); }

double SphereGrid::Phi(size_t j) const
{ return 2.0 * M_PI * j / (double)nPh; }

double SphereGrid::i(double theta) const
// { return nTh * theta / M_PI - 0.5; }
{ return (nTh - 1.0) * theta / M_PI; }

double SphereGrid::j(double phi) const
{ return nPh * phi / (2.0 * M_PI); }

size_t SphereGrid::Index(size_t i, size_t j)
{ return i + j * nTh; }

Tensor3 SphereGrid::C(size_t i, size_t j)
{
    double theta = Theta(i);
    double phi = Phi(j);
    return Tensor3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}