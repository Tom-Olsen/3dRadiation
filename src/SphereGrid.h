#ifndef __INCLUDE_GUARD_SphereGrid_h__
#define __INCLUDE_GUARD_SphereGrid_h__
#include <math.h>
#include "TensorTypes.hh"



struct SphereGrid
{
    size_t nTh;
    size_t nPh;
    size_t nDir;

    SphereGrid();
    SphereGrid(size_t nTh, size_t nPh);
    double Theta(size_t i) const;
    double Phi(size_t j) const;
    double i(double theta) const;
    double j(double phi) const;
    size_t Index(size_t i, size_t j);
    Tensor3 C(size_t i, size_t j);
};
#endif //__INCLUDE_GUARD_SphereGrid_h__