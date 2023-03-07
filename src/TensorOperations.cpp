#include "TensorOperations.h"



double Norm2(const Tensor3& vector, const Tensor3x3& gamma_ll)
{
    double norm2 = 0;
    for(int i=1; i<4; i++)
    for(int j=1; j<4; j++)
        norm2 += gamma_ll[{i,j}] * vector[i] * vector[j];
    return norm2;
}



double Norm2(const Tensor4& vector, const Tensor4x4& g_ll)
{
    double norm2 = 0;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
        norm2 += g_ll[{i,j}] * vector[i] * vector[j];
    return norm2;
}



Tensor4 NullNormalize(const Tensor4& vector, const Tensor4x4& g_ll)
{
    double a = 0;
    for(int i=1; i<4; i++)
    for(int j=1; j<4; j++)
        a += g_ll[rank2Indices{i,j}] * vector[i] * vector[j];
    double b = 0;
    for(int i=1; i<4; i++)
        b += g_ll[rank2Indices{0,i}] * vector[0] * vector[i];
    double c = g_ll[rank2Indices{0,0}] * vector[0] * vector[0];
    double d = -b / a + sqrt((b * b) / (a * a) - c / a );

    return Tensor4(vector[0], vector[1] * d, vector[2] * d, vector[3] * d);
}



Tensor3 TransformIFtoLF(const Tensor3& vector, const Tensor4x4& tetrad)
{
    return Tensor3
    (tetrad[{1,1}] * vector[1] + tetrad[{1,2}] * vector[2] + tetrad[{1,3}] * vector[3],
     tetrad[{2,1}] * vector[1] + tetrad[{2,2}] * vector[2] + tetrad[{2,3}] * vector[3],
     tetrad[{3,1}] * vector[1] + tetrad[{3,2}] * vector[2] + tetrad[{3,3}] * vector[3]);
}
Tensor3 TransformLFtoIF(const Tensor3& vector, const Tensor4x4& tetradInverse)
{
    return Tensor3
    (tetradInverse[{1,1}] * vector[1] + tetradInverse[{1,2}] * vector[2] + tetradInverse[{1,3}] * vector[3],
     tetradInverse[{2,1}] * vector[1] + tetradInverse[{2,2}] * vector[2] + tetradInverse[{2,3}] * vector[3],
     tetradInverse[{3,1}] * vector[1] + tetradInverse[{3,2}] * vector[2] + tetradInverse[{3,3}] * vector[3]);
}



Tensor4x4 TransformIFtoLF(const Tensor4x4& tensor, const Tensor4x4& tetrad)
{
    Tensor4x4 result(0.0);
    for(int a=0; a<4; a++)
        for(int b=0; b<4; b++)
            for(int A=0; A<4; A++)
                for(int B=0; B<4; B++)
                    result[{a,b}] += tensor[{A,B}] * tetrad[{a,A}] * tetrad[{b,B}];
    return result;
}
Tensor4x4 TransformLFtoIF(const Tensor4x4& tensor, const Tensor4x4& tetradInverse)
{
    Tensor4x4 result(0.0);
    for(int a=0; a<4; a++)
        for(int b=0; b<4; b++)
            for(int A=0; A<4; A++)
                for(int B=0; B<4; B++)
                    result[{a,b}] += tensor[{A,B}] * tetradInverse[{a,A}] * tetradInverse[{b,B}];
    return result;
}



template<class FrameIn, class FrameOut>
Tensor3 Vec3ObservedByEulObs(const Tensor4& u, const Coord& xyz, Metric& metric)
{
    if constexpr(std::is_same<FrameIn,IF>::value && std::is_same<FrameOut,IF>::value)
    {// IF -> IF
        double alpha = metric.GetAlpha(xyz);
        return Tensor3(u[1] / alpha, u[2] / alpha, u[3] / alpha);
    }
    if constexpr(std::is_same<FrameIn,IF>::value && std::is_same<FrameOut,LF>::value)
    {// IF -> LF
        double alpha = metric.GetAlpha(xyz);
        Tensor3 v(u[1] / alpha, u[2] / alpha, u[3] / alpha);
        return TransformIFtoLF(v, metric.GetTetrad(xyz));
    }
    if constexpr(std::is_same<FrameIn,LF>::value && std::is_same<FrameOut,IF>::value)
    {// LF -> IF
        double alpha = metric.GetAlpha(xyz);
        Tensor3 beta_u = metric.GetBeta_u(xyz);
        Tensor3 v((u[1] + beta_u[1]) / alpha, (u[2] + beta_u[2]) / alpha, (u[3] + beta_u[3]) / alpha);
        return TransformLFtoIF(v, metric.GetTetradInverse(xyz));
    }
    if constexpr(std::is_same<FrameIn,LF>::value && std::is_same<FrameOut,LF>::value)
    {// LF -> LF
        double alpha = metric.GetAlpha(xyz);
        Tensor3 beta_u = metric.GetBeta_u(xyz);
        return Tensor3((u[1] + beta_u[1]) / alpha, (u[2] + beta_u[2]) / alpha, (u[3] + beta_u[3]) / alpha);
    }
}



template Tensor3 Vec3ObservedByEulObs<IF,IF>(const Tensor4& u, const Coord& xyz, Metric& metric);
template Tensor3 Vec3ObservedByEulObs<IF,LF>(const Tensor4& u, const Coord& xyz, Metric& metric);
template Tensor3 Vec3ObservedByEulObs<LF,IF>(const Tensor4& u, const Coord& xyz, Metric& metric);
template Tensor3 Vec3ObservedByEulObs<LF,LF>(const Tensor4& u, const Coord& xyz, Metric& metric);