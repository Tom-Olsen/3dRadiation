#include <iostream>
#include "../src/Spacetimes.h"
#include "../src/TensorOperations.h"
using namespace std;


void Tetrad()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(1,1,1);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);

    Coord xyz(1.5, 2, 2.5);
    Tensor4x4 g_ll = metric.GetMetric_ll(xyz);
    Tensor4x4 tetrad = metric.GetTetrad(xyz);
    Tensor4x4 eta(0);

    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
    for(int I=0; I<4; I++)
    for(int J=0; J<4; J++)
        eta[{i,j}] += tetrad[{I,i}] * tetrad[{J,j}] * g_ll[{I,J}];

    g_ll.Print("  g_ll",1);
    tetrad.Print("tetrad",1);
    eta.Print("   eta",1);
}



void PhotonVelocity()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(1,1,1);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.0);
    
    Coord xyz(1.5, 2, 2.5);
    double alpha = metric.GetAlpha(xyz);
    Tensor4x4 g_ll = metric.GetMetric_ll(xyz);
    Tensor3x3 eta3x3_ll = metric.GetMinkowskiGamma_ll(xyz);
    Tensor4x4 eta4x4_ll = metric.GetMinkowskiMetric_ll(xyz);
    Tensor3x3 gamma_ll = metric.GetGamma_ll(xyz);

    // Four velocity defined in LF and transformed to three velocity in IF:
    {
        Tensor4 uLF(1,0.2,0.4,0.6);
        uLF = NullNormalize(uLF,g_ll);
        Tensor3 vIF = Vec3ObservedByEulObs<LF,IF>(uLF,xyz,metric);

        uLF.Print(" uLF ");
        PrintDouble(Norm2(uLF,g_ll),"|uLF|");
        vIF.Print(" vIF ");
        PrintDouble(Norm2(vIF,eta3x3_ll),"|vIF|");
    }cout << endl;
    // Four velocity defined in IF and transformed to three velocity in LF:
    {
        Tensor4 uIF(1*alpha,0.2*alpha,0.4*alpha,0.6*alpha);
        Tensor3 vLF = Vec3ObservedByEulObs<IF,LF>(uIF, xyz, metric);

        uIF.Print(" uIF ");
        PrintDouble(Norm2(uIF,eta4x4_ll),"|uIF|");
        vLF.Print(" vLF ");
        PrintDouble(Norm2(vLF,gamma_ll),"|vLF|");
    }cout << endl;
    // Four velocity defined in LF and transformed to three velocity in LF:
    {
        Tensor4 uLF(1,0.2,0.4,0.6);
        uLF = NullNormalize(uLF,g_ll);
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF, xyz, metric);

        uLF.Print(" uLF ");
        PrintDouble(Norm2(uLF,g_ll),"|uLF|");
        vLF.Print(" vLF ");
        PrintDouble(Norm2(vLF,gamma_ll),"|vLF|");
    }
}



int main()
{
    Tetrad();
    PhotonVelocity();
}