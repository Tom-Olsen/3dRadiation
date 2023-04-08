#ifndef __INCLUDE_GUARD_TensorOperations_h__
#define __INCLUDE_GUARD_TensorOperations_h__
#include "TensorTypes.hh"   // General relativity tensors.
#include "Spacetimes.h"     // Metric data.

double Norm2(const Tensor3& vector, const Tensor3x3& gamma_ll);
double Norm2(const Tensor4& vector, const Tensor4x4& g_ll);
Tensor4 NullNormalize(const Tensor4& vector, const Tensor4x4& g_ll);

Tensor3 TransformIFtoLF(const Tensor3& vector, const Tensor4x4& tetrad);
Tensor3 TransformLFtoIF(const Tensor3& vector, const Tensor4x4& tetradInverse);

Tensor4x4 TransformIFtoLF(const Tensor4x4& tensor, const Tensor4x4& tetrad);
Tensor4x4 TransformLFtoIF(const Tensor4x4& tensor, const Tensor4x4& tetradInverse);

template<class FrameIn, class FrameOut>
Tensor3 Vec3ObservedByEulObs(const Tensor4& u, const Coord& xyz, Metric& metric);

#endif //__INCLUDE_GUARD_TensorOperations_h__