#ifndef __INCLUDE_GUARD_Metric_h__
#define __INCLUDE_GUARD_Metric_h__
#include <string>
#include "ControlFlow.hh"       // Template arguments and profiling macros.
#include "Utility.hh"           // Utility functions.
#include "TensorTypes.hh"       // General relativity tensors.
#include "Grid.h"               // Numerical Grid and mapping to physical domain.
#include "Interpolation.hh"     // Basic interpolation schemes.



class Metric
{
public:
    // Grid Data:
    Grid& grid;
    double m = 1.0;
    double a = 0.0;

// protected:
public:
    // Metrik Data:
    RealBuffer g00_ll;    RealBuffer g00_uu;
    RealBuffer g01_ll;    RealBuffer g01_uu;
    RealBuffer g02_ll;    RealBuffer g02_uu;
    RealBuffer g03_ll;    RealBuffer g03_uu;
    RealBuffer g11_ll;    RealBuffer g11_uu;
    RealBuffer g12_ll;    RealBuffer g12_uu;
    RealBuffer g13_ll;    RealBuffer g13_uu;
    RealBuffer g22_ll;    RealBuffer g22_uu;
    RealBuffer g23_ll;    RealBuffer g23_uu;
    RealBuffer g33_ll;    RealBuffer g33_uu;
    RealBuffer d0_g00_lll;    RealBuffer d1_g00_lll;    RealBuffer d2_g00_lll;    RealBuffer d3_g00_lll;
    RealBuffer d0_g01_lll;    RealBuffer d1_g01_lll;    RealBuffer d2_g01_lll;    RealBuffer d3_g01_lll;
    RealBuffer d0_g02_lll;    RealBuffer d1_g02_lll;    RealBuffer d2_g02_lll;    RealBuffer d3_g02_lll;
    RealBuffer d0_g03_lll;    RealBuffer d1_g03_lll;    RealBuffer d2_g03_lll;    RealBuffer d3_g03_lll;
    RealBuffer d0_g11_lll;    RealBuffer d1_g11_lll;    RealBuffer d2_g11_lll;    RealBuffer d3_g11_lll;
    RealBuffer d0_g12_lll;    RealBuffer d1_g12_lll;    RealBuffer d2_g12_lll;    RealBuffer d3_g12_lll;
    RealBuffer d0_g13_lll;    RealBuffer d1_g13_lll;    RealBuffer d2_g13_lll;    RealBuffer d3_g13_lll;
    RealBuffer d0_g22_lll;    RealBuffer d1_g22_lll;    RealBuffer d2_g22_lll;    RealBuffer d3_g22_lll;
    RealBuffer d0_g23_lll;    RealBuffer d1_g23_lll;    RealBuffer d2_g23_lll;    RealBuffer d3_g23_lll;
    RealBuffer d0_g33_lll;    RealBuffer d1_g33_lll;    RealBuffer d2_g33_lll;    RealBuffer d3_g33_lll;
    RealBuffer d0_g00_luu;    RealBuffer d1_g00_luu;    RealBuffer d2_g00_luu;    RealBuffer d3_g00_luu;
    RealBuffer d0_g01_luu;    RealBuffer d1_g01_luu;    RealBuffer d2_g01_luu;    RealBuffer d3_g01_luu;
    RealBuffer d0_g02_luu;    RealBuffer d1_g02_luu;    RealBuffer d2_g02_luu;    RealBuffer d3_g02_luu;
    RealBuffer d0_g03_luu;    RealBuffer d1_g03_luu;    RealBuffer d2_g03_luu;    RealBuffer d3_g03_luu;
    RealBuffer d0_g11_luu;    RealBuffer d1_g11_luu;    RealBuffer d2_g11_luu;    RealBuffer d3_g11_luu;
    RealBuffer d0_g12_luu;    RealBuffer d1_g12_luu;    RealBuffer d2_g12_luu;    RealBuffer d3_g12_luu;
    RealBuffer d0_g13_luu;    RealBuffer d1_g13_luu;    RealBuffer d2_g13_luu;    RealBuffer d3_g13_luu;
    RealBuffer d0_g22_luu;    RealBuffer d1_g22_luu;    RealBuffer d2_g22_luu;    RealBuffer d3_g22_luu;
    RealBuffer d0_g23_luu;    RealBuffer d1_g23_luu;    RealBuffer d2_g23_luu;    RealBuffer d3_g23_luu;
    RealBuffer d0_g33_luu;    RealBuffer d1_g33_luu;    RealBuffer d2_g33_luu;    RealBuffer d3_g33_luu;
    // ADM Data:
    RealBuffer alpha;
    RealBuffer beta1_u;    RealBuffer beta1_l;
    RealBuffer beta2_u;    RealBuffer beta2_l;
    RealBuffer beta3_u;    RealBuffer beta3_l;
    RealBuffer gamma11_ll;    RealBuffer gamma11_uu;
    RealBuffer gamma12_ll;    RealBuffer gamma12_uu;
    RealBuffer gamma13_ll;    RealBuffer gamma13_uu;
    RealBuffer gamma22_ll;    RealBuffer gamma22_uu;
    RealBuffer gamma23_ll;    RealBuffer gamma23_uu;
    RealBuffer gamma33_ll;    RealBuffer gamma33_uu;
    RealBuffer d1_alpha_l;
    RealBuffer d2_alpha_l;
    RealBuffer d3_alpha_l;
    RealBuffer d1_beta1_lu;    RealBuffer d1_beta2_lu;    RealBuffer d1_beta3_lu;
    RealBuffer d2_beta1_lu;    RealBuffer d2_beta2_lu;    RealBuffer d2_beta3_lu;
    RealBuffer d3_beta1_lu;    RealBuffer d3_beta2_lu;    RealBuffer d3_beta3_lu;
    RealBuffer d1_beta1_ll;    RealBuffer d1_beta2_ll;    RealBuffer d1_beta3_ll;
    RealBuffer d2_beta1_ll;    RealBuffer d2_beta2_ll;    RealBuffer d2_beta3_ll;
    RealBuffer d3_beta1_ll;    RealBuffer d3_beta2_ll;    RealBuffer d3_beta3_ll;
    RealBuffer d1_gamma11_lll;    RealBuffer d2_gamma11_lll;    RealBuffer d3_gamma11_lll;
    RealBuffer d1_gamma12_lll;    RealBuffer d2_gamma12_lll;    RealBuffer d3_gamma12_lll;
    RealBuffer d1_gamma13_lll;    RealBuffer d2_gamma13_lll;    RealBuffer d3_gamma13_lll;
    RealBuffer d1_gamma22_lll;    RealBuffer d2_gamma22_lll;    RealBuffer d3_gamma22_lll;
    RealBuffer d1_gamma23_lll;    RealBuffer d2_gamma23_lll;    RealBuffer d3_gamma23_lll;
    RealBuffer d1_gamma33_lll;    RealBuffer d2_gamma33_lll;    RealBuffer d3_gamma33_lll;
    RealBuffer d1_gamma11_luu;    RealBuffer d2_gamma11_luu;    RealBuffer d3_gamma11_luu;
    RealBuffer d1_gamma12_luu;    RealBuffer d2_gamma12_luu;    RealBuffer d3_gamma12_luu;
    RealBuffer d1_gamma13_luu;    RealBuffer d2_gamma13_luu;    RealBuffer d3_gamma13_luu;
    RealBuffer d1_gamma22_luu;    RealBuffer d2_gamma22_luu;    RealBuffer d3_gamma22_luu;
    RealBuffer d1_gamma23_luu;    RealBuffer d2_gamma23_luu;    RealBuffer d3_gamma23_luu;
    RealBuffer d1_gamma33_luu;    RealBuffer d2_gamma33_luu;    RealBuffer d3_gamma33_luu;
    RealBuffer K11_ll;
    RealBuffer K12_ll;
    RealBuffer K13_ll;
    RealBuffer K22_ll;
    RealBuffer K23_ll;
    RealBuffer K33_ll;
    // Tetrad Data:
    RealBuffer tetrad00_ul;    RealBuffer tetrad01_ul;    RealBuffer tetrad02_ul;    RealBuffer tetrad03_ul;
    RealBuffer tetrad10_ul;    RealBuffer tetrad11_ul;    RealBuffer tetrad12_ul;    RealBuffer tetrad13_ul;
    RealBuffer tetrad20_ul;    RealBuffer tetrad21_ul;    RealBuffer tetrad22_ul;    RealBuffer tetrad23_ul;
    RealBuffer tetrad30_ul;    RealBuffer tetrad31_ul;    RealBuffer tetrad32_ul;    RealBuffer tetrad33_ul;
    RealBuffer tetrad00_lu;    RealBuffer tetrad01_lu;    RealBuffer tetrad02_lu;    RealBuffer tetrad03_lu;
    RealBuffer tetrad10_lu;    RealBuffer tetrad11_lu;    RealBuffer tetrad12_lu;    RealBuffer tetrad13_lu;
    RealBuffer tetrad20_lu;    RealBuffer tetrad21_lu;    RealBuffer tetrad22_lu;    RealBuffer tetrad23_lu;
    RealBuffer tetrad30_lu;    RealBuffer tetrad31_lu;    RealBuffer tetrad32_lu;    RealBuffer tetrad33_lu;

public:
    // Constructors/Destructor:
    Metric(Grid& grid_, double m_, double a_);

    virtual std::string Name();

    // Initialization:
    virtual Tensor4x4 MetricFunction(const Coord& xyz);
    void InitializeMetricOnGrid();
    void InitializeBoostedTetradOnGrid();
    template<int k>
    Tensor4x4 MetricDeriv(const Coord& xyz);
    template<int k>
    Tensor4x4 InverseMetricDeriv(const Coord& xyz);
    void InitializeMetricDerivativesOnGrid();
    void InitializeAdmComponentsOnGrid();
    double InterpolateArrayTo_ijk(const RealBuffer& array, const Coord& ijk);
    double InterpolateArrayTo_ijk(const RealBuffer& array, double i, double j, double k);

public:
    // Boolean checks:
    virtual bool InsideBH(const Coord& xyz);

    // Tensor getters:
    Tensor4 uEulObs(int ijk);
    Tensor4 uEulObs(const Coord& xyz);
    Tensor4x4 GetMetric_ll(int ijk);
    Tensor4x4 GetMetric_ll(const Coord& xyz);
    Tensor4x4 GetMetric_uu(int ijk);
    Tensor4x4 GetMetric_uu(const Coord& xyz);
    Tensor4x4 GetMinkowskiMetric_ll(int ijk);
    Tensor4x4 GetMinkowskiMetric_ll(const Coord& xyz);
    Tensor4x4 GetMinkowskiMetric_uu(int ijk);
    Tensor4x4 GetMinkowskiMetric_uu(const Coord& xyz);
    Tensor4x4x4 GetDerivMetric_lll(int ijk);
    Tensor4x4x4 GetDerivMetric_lll(const Coord& xyz);
    Tensor4x4x4 GetDerivMetric_luu(int ijk);
    Tensor4x4x4 GetDerivMetric_luu(const Coord& xyz);
    Tensor4x4 GetTetrad(int ijk);
    Tensor4x4 GetTetrad(const Coord& xyz);
    Tensor4x4 GetTetradInverse(int ijk);
    Tensor4x4 GetTetradInverse(const Coord& xyz);

    // ADM getters:
    double GetAlpha(int ijk);
    double GetAlpha(const Coord& xyz);
    Tensor3 GetBeta_u(int ijk);
    Tensor3 GetBeta_u(const Coord& xyz);
    Tensor3 GetBeta_l(int ijk);
    Tensor3 GetBeta_l(const Coord& xyz);
    Tensor3x3 GetGamma_ll(int ijk);
    Tensor3x3 GetGamma_ll(const Coord& xyz);
    Tensor3x3 GetGamma_uu(int ijk);
    Tensor3x3 GetGamma_uu(const Coord& xyz);
    Tensor3x3 GetMinkowskiGamma_ll(int ijk);
    Tensor3x3 GetMinkowskiGamma_ll(const Coord& xyz);
    Tensor3x3 GetMinkowskiGamma_uu(int ijk);
    Tensor3x3 GetMinkowskiGamma_uu(const Coord& xyz);

    Tensor3 GetDerivAlpha_l(int ijk);
    Tensor3 GetDerivAlpha_l(const Coord& xyz);
    Tensor3x3 GetDerivBeta_lu(int ijk);
    Tensor3x3 GetDerivBeta_lu(const Coord& xyz);
    Tensor3x3 GetDerivBeta_ll(int ijk);
    Tensor3x3 GetDerivBeta_ll(const Coord& xyz);
    Tensor3x3x3 GetDerivGamma_lll(int ijk);
    Tensor3x3x3 GetDerivGamma_lll(const Coord& xyz);
    Tensor3x3x3 GetDerivGamma_luu(int ijk);
    Tensor3x3x3 GetDerivGamma_luu(const Coord& xyz);
    Tensor3x3 GetExtrCurv_ll(int ijk);
    Tensor3x3 GetExtrCurv_ll(const Coord& xyz);
};
#endif //__INCLUDE_GUARD_Metric_h__