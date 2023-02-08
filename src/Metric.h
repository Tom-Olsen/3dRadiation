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
    double* g00_ll;    double* g00_uu;
    double* g01_ll;    double* g01_uu;
    double* g02_ll;    double* g02_uu;
    double* g03_ll;    double* g03_uu;
    double* g11_ll;    double* g11_uu;
    double* g12_ll;    double* g12_uu;
    double* g13_ll;    double* g13_uu;
    double* g22_ll;    double* g22_uu;
    double* g23_ll;    double* g23_uu;
    double* g33_ll;    double* g33_uu;
    double* d0_g00_lll;    double* d1_g00_lll;    double* d2_g00_lll;    double* d3_g00_lll;
    double* d0_g01_lll;    double* d1_g01_lll;    double* d2_g01_lll;    double* d3_g01_lll;
    double* d0_g02_lll;    double* d1_g02_lll;    double* d2_g02_lll;    double* d3_g02_lll;
    double* d0_g03_lll;    double* d1_g03_lll;    double* d2_g03_lll;    double* d3_g03_lll;
    double* d0_g11_lll;    double* d1_g11_lll;    double* d2_g11_lll;    double* d3_g11_lll;
    double* d0_g12_lll;    double* d1_g12_lll;    double* d2_g12_lll;    double* d3_g12_lll;
    double* d0_g13_lll;    double* d1_g13_lll;    double* d2_g13_lll;    double* d3_g13_lll;
    double* d0_g22_lll;    double* d1_g22_lll;    double* d2_g22_lll;    double* d3_g22_lll;
    double* d0_g23_lll;    double* d1_g23_lll;    double* d2_g23_lll;    double* d3_g23_lll;
    double* d0_g33_lll;    double* d1_g33_lll;    double* d2_g33_lll;    double* d3_g33_lll;
    double* d0_g00_luu;    double* d1_g00_luu;    double* d2_g00_luu;    double* d3_g00_luu;
    double* d0_g01_luu;    double* d1_g01_luu;    double* d2_g01_luu;    double* d3_g01_luu;
    double* d0_g02_luu;    double* d1_g02_luu;    double* d2_g02_luu;    double* d3_g02_luu;
    double* d0_g03_luu;    double* d1_g03_luu;    double* d2_g03_luu;    double* d3_g03_luu;
    double* d0_g11_luu;    double* d1_g11_luu;    double* d2_g11_luu;    double* d3_g11_luu;
    double* d0_g12_luu;    double* d1_g12_luu;    double* d2_g12_luu;    double* d3_g12_luu;
    double* d0_g13_luu;    double* d1_g13_luu;    double* d2_g13_luu;    double* d3_g13_luu;
    double* d0_g22_luu;    double* d1_g22_luu;    double* d2_g22_luu;    double* d3_g22_luu;
    double* d0_g23_luu;    double* d1_g23_luu;    double* d2_g23_luu;    double* d3_g23_luu;
    double* d0_g33_luu;    double* d1_g33_luu;    double* d2_g33_luu;    double* d3_g33_luu;
    // ADM Data:
    double* alpha;
    double* beta1_u;    double* beta1_l;
    double* beta2_u;    double* beta2_l;
    double* beta3_u;    double* beta3_l;
    double* gamma11_ll;    double* gamma11_uu;
    double* gamma12_ll;    double* gamma12_uu;
    double* gamma13_ll;    double* gamma13_uu;
    double* gamma22_ll;    double* gamma22_uu;
    double* gamma23_ll;    double* gamma23_uu;
    double* gamma33_ll;    double* gamma33_uu;
    double* d1_alpha_l;
    double* d2_alpha_l;
    double* d3_alpha_l;
    double* d1_beta1_lu;    double* d1_beta2_lu;    double* d1_beta3_lu;
    double* d2_beta1_lu;    double* d2_beta2_lu;    double* d2_beta3_lu;
    double* d3_beta1_lu;    double* d3_beta2_lu;    double* d3_beta3_lu;
    double* d1_beta1_ll;    double* d1_beta2_ll;    double* d1_beta3_ll;
    double* d2_beta1_ll;    double* d2_beta2_ll;    double* d2_beta3_ll;
    double* d3_beta1_ll;    double* d3_beta2_ll;    double* d3_beta3_ll;
    double* d1_gamma11_lll;    double* d2_gamma11_lll;    double* d3_gamma11_lll;
    double* d1_gamma12_lll;    double* d2_gamma12_lll;    double* d3_gamma12_lll;
    double* d1_gamma13_lll;    double* d2_gamma13_lll;    double* d3_gamma13_lll;
    double* d1_gamma22_lll;    double* d2_gamma22_lll;    double* d3_gamma22_lll;
    double* d1_gamma23_lll;    double* d2_gamma23_lll;    double* d3_gamma23_lll;
    double* d1_gamma33_lll;    double* d2_gamma33_lll;    double* d3_gamma33_lll;
    double* d1_gamma11_luu;    double* d2_gamma11_luu;    double* d3_gamma11_luu;
    double* d1_gamma12_luu;    double* d2_gamma12_luu;    double* d3_gamma12_luu;
    double* d1_gamma13_luu;    double* d2_gamma13_luu;    double* d3_gamma13_luu;
    double* d1_gamma22_luu;    double* d2_gamma22_luu;    double* d3_gamma22_luu;
    double* d1_gamma23_luu;    double* d2_gamma23_luu;    double* d3_gamma23_luu;
    double* d1_gamma33_luu;    double* d2_gamma33_luu;    double* d3_gamma33_luu;
    double* K11_ll;
    double* K12_ll;
    double* K13_ll;
    double* K22_ll;
    double* K23_ll;
    double* K33_ll;
    // Tetrad Data:
    double* tetrad00_ul;    double* tetrad01_ul;    double* tetrad02_ul;    double* tetrad03_ul;
    double* tetrad10_ul;    double* tetrad11_ul;    double* tetrad12_ul;    double* tetrad13_ul;
    double* tetrad20_ul;    double* tetrad21_ul;    double* tetrad22_ul;    double* tetrad23_ul;
    double* tetrad30_ul;    double* tetrad31_ul;    double* tetrad32_ul;    double* tetrad33_ul;

public:
    // Constructors/Destructor:
    Metric(Grid& grid_, double m_, double a_);
    ~Metric();

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
    double InterpolateArrayTo_ijk(double* array, const Coord& ijk);
    double InterpolateArrayTo_ijk(double* array, double i, double j, double k);

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