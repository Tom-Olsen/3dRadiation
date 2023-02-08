#include "Metric.h"



// Constructor/Desctructor
Metric::Metric(Grid& grid_, double m_, double a_) : grid(grid_), m(m_), a(a_)
{
    // Metrik Data:
    g00_ll = new double[grid.nxyz]();    g00_uu = new double[grid.nxyz]();
    g01_ll = new double[grid.nxyz]();    g01_uu = new double[grid.nxyz]();
    g02_ll = new double[grid.nxyz]();    g02_uu = new double[grid.nxyz]();
    g03_ll = new double[grid.nxyz]();    g03_uu = new double[grid.nxyz]();
    g11_ll = new double[grid.nxyz]();    g11_uu = new double[grid.nxyz]();
    g12_ll = new double[grid.nxyz]();    g12_uu = new double[grid.nxyz]();
    g13_ll = new double[grid.nxyz]();    g13_uu = new double[grid.nxyz]();
    g22_ll = new double[grid.nxyz]();    g22_uu = new double[grid.nxyz]();
    g23_ll = new double[grid.nxyz]();    g23_uu = new double[grid.nxyz]();
    g33_ll = new double[grid.nxyz]();    g33_uu = new double[grid.nxyz]();
    d0_g00_lll = new double[grid.nxyz]();    d1_g00_lll = new double[grid.nxyz]();    d2_g00_lll = new double[grid.nxyz]();    d3_g00_lll = new double[grid.nxyz]();
    d0_g01_lll = new double[grid.nxyz]();    d1_g01_lll = new double[grid.nxyz]();    d2_g01_lll = new double[grid.nxyz]();    d3_g01_lll = new double[grid.nxyz]();
    d0_g02_lll = new double[grid.nxyz]();    d1_g02_lll = new double[grid.nxyz]();    d2_g02_lll = new double[grid.nxyz]();    d3_g02_lll = new double[grid.nxyz]();
    d0_g03_lll = new double[grid.nxyz]();    d1_g03_lll = new double[grid.nxyz]();    d2_g03_lll = new double[grid.nxyz]();    d3_g03_lll = new double[grid.nxyz]();
    d0_g11_lll = new double[grid.nxyz]();    d1_g11_lll = new double[grid.nxyz]();    d2_g11_lll = new double[grid.nxyz]();    d3_g11_lll = new double[grid.nxyz]();
    d0_g12_lll = new double[grid.nxyz]();    d1_g12_lll = new double[grid.nxyz]();    d2_g12_lll = new double[grid.nxyz]();    d3_g12_lll = new double[grid.nxyz]();
    d0_g13_lll = new double[grid.nxyz]();    d1_g13_lll = new double[grid.nxyz]();    d2_g13_lll = new double[grid.nxyz]();    d3_g13_lll = new double[grid.nxyz]();
    d0_g22_lll = new double[grid.nxyz]();    d1_g22_lll = new double[grid.nxyz]();    d2_g22_lll = new double[grid.nxyz]();    d3_g22_lll = new double[grid.nxyz]();
    d0_g23_lll = new double[grid.nxyz]();    d1_g23_lll = new double[grid.nxyz]();    d2_g23_lll = new double[grid.nxyz]();    d3_g23_lll = new double[grid.nxyz]();
    d0_g33_lll = new double[grid.nxyz]();    d1_g33_lll = new double[grid.nxyz]();    d2_g33_lll = new double[grid.nxyz]();    d3_g33_lll = new double[grid.nxyz]();
    d0_g00_luu = new double[grid.nxyz]();    d1_g00_luu = new double[grid.nxyz]();    d2_g00_luu = new double[grid.nxyz]();    d3_g00_luu = new double[grid.nxyz]();
    d0_g01_luu = new double[grid.nxyz]();    d1_g01_luu = new double[grid.nxyz]();    d2_g01_luu = new double[grid.nxyz]();    d3_g01_luu = new double[grid.nxyz]();
    d0_g02_luu = new double[grid.nxyz]();    d1_g02_luu = new double[grid.nxyz]();    d2_g02_luu = new double[grid.nxyz]();    d3_g02_luu = new double[grid.nxyz]();
    d0_g03_luu = new double[grid.nxyz]();    d1_g03_luu = new double[grid.nxyz]();    d2_g03_luu = new double[grid.nxyz]();    d3_g03_luu = new double[grid.nxyz]();
    d0_g11_luu = new double[grid.nxyz]();    d1_g11_luu = new double[grid.nxyz]();    d2_g11_luu = new double[grid.nxyz]();    d3_g11_luu = new double[grid.nxyz]();
    d0_g12_luu = new double[grid.nxyz]();    d1_g12_luu = new double[grid.nxyz]();    d2_g12_luu = new double[grid.nxyz]();    d3_g12_luu = new double[grid.nxyz]();
    d0_g13_luu = new double[grid.nxyz]();    d1_g13_luu = new double[grid.nxyz]();    d2_g13_luu = new double[grid.nxyz]();    d3_g13_luu = new double[grid.nxyz]();
    d0_g22_luu = new double[grid.nxyz]();    d1_g22_luu = new double[grid.nxyz]();    d2_g22_luu = new double[grid.nxyz]();    d3_g22_luu = new double[grid.nxyz]();
    d0_g23_luu = new double[grid.nxyz]();    d1_g23_luu = new double[grid.nxyz]();    d2_g23_luu = new double[grid.nxyz]();    d3_g23_luu = new double[grid.nxyz]();
    d0_g33_luu = new double[grid.nxyz]();    d1_g33_luu = new double[grid.nxyz]();    d2_g33_luu = new double[grid.nxyz]();    d3_g33_luu = new double[grid.nxyz]();
    // ADM Data:
    alpha = new double[grid.nxyz]();
    beta1_u = new double[grid.nxyz]();    beta1_l = new double[grid.nxyz]();
    beta2_u = new double[grid.nxyz]();    beta2_l = new double[grid.nxyz]();
    beta3_u = new double[grid.nxyz]();    beta3_l = new double[grid.nxyz]();
    gamma11_ll = new double[grid.nxyz]();    gamma11_uu = new double[grid.nxyz]();
    gamma12_ll = new double[grid.nxyz]();    gamma12_uu = new double[grid.nxyz]();
    gamma13_ll = new double[grid.nxyz]();    gamma13_uu = new double[grid.nxyz]();
    gamma22_ll = new double[grid.nxyz]();    gamma22_uu = new double[grid.nxyz]();
    gamma23_ll = new double[grid.nxyz]();    gamma23_uu = new double[grid.nxyz]();
    gamma33_ll = new double[grid.nxyz]();    gamma33_uu = new double[grid.nxyz]();
    d1_alpha_l = new double[grid.nxyz]();
    d2_alpha_l = new double[grid.nxyz]();
    d3_alpha_l = new double[grid.nxyz]();
    d1_beta1_lu = new double[grid.nxyz]();    d1_beta2_lu = new double[grid.nxyz]();    d1_beta3_lu = new double[grid.nxyz]();
    d2_beta1_lu = new double[grid.nxyz]();    d2_beta2_lu = new double[grid.nxyz]();    d2_beta3_lu = new double[grid.nxyz]();
    d3_beta1_lu = new double[grid.nxyz]();    d3_beta2_lu = new double[grid.nxyz]();    d3_beta3_lu = new double[grid.nxyz]();
    d1_beta1_ll = new double[grid.nxyz]();    d1_beta2_ll = new double[grid.nxyz]();    d1_beta3_ll = new double[grid.nxyz]();
    d2_beta1_ll = new double[grid.nxyz]();    d2_beta2_ll = new double[grid.nxyz]();    d2_beta3_ll = new double[grid.nxyz]();
    d3_beta1_ll = new double[grid.nxyz]();    d3_beta2_ll = new double[grid.nxyz]();    d3_beta3_ll = new double[grid.nxyz]();
    d1_gamma11_lll = new double[grid.nxyz]();    d2_gamma11_lll = new double[grid.nxyz]();    d3_gamma11_lll = new double[grid.nxyz]();
    d1_gamma12_lll = new double[grid.nxyz]();    d2_gamma12_lll = new double[grid.nxyz]();    d3_gamma12_lll = new double[grid.nxyz]();
    d1_gamma13_lll = new double[grid.nxyz]();    d2_gamma13_lll = new double[grid.nxyz]();    d3_gamma13_lll = new double[grid.nxyz]();
    d1_gamma22_lll = new double[grid.nxyz]();    d2_gamma22_lll = new double[grid.nxyz]();    d3_gamma22_lll = new double[grid.nxyz]();
    d1_gamma23_lll = new double[grid.nxyz]();    d2_gamma23_lll = new double[grid.nxyz]();    d3_gamma23_lll = new double[grid.nxyz]();
    d1_gamma33_lll = new double[grid.nxyz]();    d2_gamma33_lll = new double[grid.nxyz]();    d3_gamma33_lll = new double[grid.nxyz]();
    d1_gamma11_luu = new double[grid.nxyz]();    d2_gamma11_luu = new double[grid.nxyz]();    d3_gamma11_luu = new double[grid.nxyz]();
    d1_gamma12_luu = new double[grid.nxyz]();    d2_gamma12_luu = new double[grid.nxyz]();    d3_gamma12_luu = new double[grid.nxyz]();
    d1_gamma13_luu = new double[grid.nxyz]();    d2_gamma13_luu = new double[grid.nxyz]();    d3_gamma13_luu = new double[grid.nxyz]();
    d1_gamma22_luu = new double[grid.nxyz]();    d2_gamma22_luu = new double[grid.nxyz]();    d3_gamma22_luu = new double[grid.nxyz]();
    d1_gamma23_luu = new double[grid.nxyz]();    d2_gamma23_luu = new double[grid.nxyz]();    d3_gamma23_luu = new double[grid.nxyz]();
    d1_gamma33_luu = new double[grid.nxyz]();    d2_gamma33_luu = new double[grid.nxyz]();    d3_gamma33_luu = new double[grid.nxyz]();
    K11_ll = new double[grid.nxyz]();
    K12_ll = new double[grid.nxyz]();
    K13_ll = new double[grid.nxyz]();
    K22_ll = new double[grid.nxyz]();
    K23_ll = new double[grid.nxyz]();
    K33_ll = new double[grid.nxyz]();
    // Tetrad Data:
    tetrad00_ul = new double[grid.nxyz]();    tetrad01_ul = new double[grid.nxyz]();    tetrad02_ul = new double[grid.nxyz]();    tetrad03_ul = new double[grid.nxyz]();
    tetrad10_ul = new double[grid.nxyz]();    tetrad11_ul = new double[grid.nxyz]();    tetrad12_ul = new double[grid.nxyz]();    tetrad13_ul = new double[grid.nxyz]();
    tetrad20_ul = new double[grid.nxyz]();    tetrad21_ul = new double[grid.nxyz]();    tetrad22_ul = new double[grid.nxyz]();    tetrad23_ul = new double[grid.nxyz]();
    tetrad30_ul = new double[grid.nxyz]();    tetrad31_ul = new double[grid.nxyz]();    tetrad32_ul = new double[grid.nxyz]();    tetrad33_ul = new double[grid.nxyz]();
}
Metric::~Metric()
{
    // Metrik Data:
    delete[] g00_ll;    delete[] g00_uu;
    delete[] g01_ll;    delete[] g01_uu;
    delete[] g02_ll;    delete[] g02_uu;
    delete[] g03_ll;    delete[] g03_uu;
    delete[] g11_ll;    delete[] g11_uu;
    delete[] g12_ll;    delete[] g12_uu;
    delete[] g13_ll;    delete[] g13_uu;
    delete[] g22_ll;    delete[] g22_uu;
    delete[] g23_ll;    delete[] g23_uu;
    delete[] g33_ll;    delete[] g33_uu;
    delete[] d0_g00_lll;    delete[] d1_g00_lll;    delete[] d2_g00_lll;    delete[] d3_g00_lll;
    delete[] d0_g01_lll;    delete[] d1_g01_lll;    delete[] d2_g01_lll;    delete[] d3_g01_lll;
    delete[] d0_g02_lll;    delete[] d1_g02_lll;    delete[] d2_g02_lll;    delete[] d3_g02_lll;
    delete[] d0_g03_lll;    delete[] d1_g03_lll;    delete[] d2_g03_lll;    delete[] d3_g03_lll;
    delete[] d0_g11_lll;    delete[] d1_g11_lll;    delete[] d2_g11_lll;    delete[] d3_g11_lll;
    delete[] d0_g12_lll;    delete[] d1_g12_lll;    delete[] d2_g12_lll;    delete[] d3_g12_lll;
    delete[] d0_g13_lll;    delete[] d1_g13_lll;    delete[] d2_g13_lll;    delete[] d3_g13_lll;
    delete[] d0_g22_lll;    delete[] d1_g22_lll;    delete[] d2_g22_lll;    delete[] d3_g22_lll;
    delete[] d0_g23_lll;    delete[] d1_g23_lll;    delete[] d2_g23_lll;    delete[] d3_g23_lll;
    delete[] d0_g33_lll;    delete[] d1_g33_lll;    delete[] d2_g33_lll;    delete[] d3_g33_lll;
    delete[] d0_g00_luu;    delete[] d1_g00_luu;    delete[] d2_g00_luu;    delete[] d3_g00_luu;
    delete[] d0_g01_luu;    delete[] d1_g01_luu;    delete[] d2_g01_luu;    delete[] d3_g01_luu;
    delete[] d0_g02_luu;    delete[] d1_g02_luu;    delete[] d2_g02_luu;    delete[] d3_g02_luu;
    delete[] d0_g03_luu;    delete[] d1_g03_luu;    delete[] d2_g03_luu;    delete[] d3_g03_luu;
    delete[] d0_g11_luu;    delete[] d1_g11_luu;    delete[] d2_g11_luu;    delete[] d3_g11_luu;
    delete[] d0_g12_luu;    delete[] d1_g12_luu;    delete[] d2_g12_luu;    delete[] d3_g12_luu;
    delete[] d0_g13_luu;    delete[] d1_g13_luu;    delete[] d2_g13_luu;    delete[] d3_g13_luu;
    delete[] d0_g22_luu;    delete[] d1_g22_luu;    delete[] d2_g22_luu;    delete[] d3_g22_luu;
    delete[] d0_g23_luu;    delete[] d1_g23_luu;    delete[] d2_g23_luu;    delete[] d3_g23_luu;
    delete[] d0_g33_luu;    delete[] d1_g33_luu;    delete[] d2_g33_luu;    delete[] d3_g33_luu;
    // ADM Data:
    delete[] alpha;
    delete[] beta1_u;    delete[] beta1_l;
    delete[] beta2_u;    delete[] beta2_l;
    delete[] beta3_u;    delete[] beta3_l;
    delete[] gamma11_ll;    delete[] gamma11_uu;
    delete[] gamma12_ll;    delete[] gamma12_uu;
    delete[] gamma13_ll;    delete[] gamma13_uu;
    delete[] gamma22_ll;    delete[] gamma22_uu;
    delete[] gamma23_ll;    delete[] gamma23_uu;
    delete[] gamma33_ll;    delete[] gamma33_uu;
    delete[] d1_alpha_l;
    delete[] d2_alpha_l;
    delete[] d3_alpha_l;
    delete[] d1_beta1_lu;    delete[] d1_beta2_lu;    delete[] d1_beta3_lu;
    delete[] d2_beta1_lu;    delete[] d2_beta2_lu;    delete[] d2_beta3_lu;
    delete[] d3_beta1_lu;    delete[] d3_beta2_lu;    delete[] d3_beta3_lu;
    delete[] d1_beta1_ll;    delete[] d1_beta2_ll;    delete[] d1_beta3_ll;
    delete[] d2_beta1_ll;    delete[] d2_beta2_ll;    delete[] d2_beta3_ll;
    delete[] d3_beta1_ll;    delete[] d3_beta2_ll;    delete[] d3_beta3_ll;
    delete[] d1_gamma11_lll;    delete[] d2_gamma11_lll;    delete[] d3_gamma11_lll;
    delete[] d1_gamma12_lll;    delete[] d2_gamma12_lll;    delete[] d3_gamma12_lll;
    delete[] d1_gamma13_lll;    delete[] d2_gamma13_lll;    delete[] d3_gamma13_lll;
    delete[] d1_gamma22_lll;    delete[] d2_gamma22_lll;    delete[] d3_gamma22_lll;
    delete[] d1_gamma23_lll;    delete[] d2_gamma23_lll;    delete[] d3_gamma23_lll;
    delete[] d1_gamma33_lll;    delete[] d2_gamma33_lll;    delete[] d3_gamma33_lll;
    delete[] d1_gamma11_luu;    delete[] d2_gamma11_luu;    delete[] d3_gamma11_luu;
    delete[] d1_gamma12_luu;    delete[] d2_gamma12_luu;    delete[] d3_gamma12_luu;
    delete[] d1_gamma13_luu;    delete[] d2_gamma13_luu;    delete[] d3_gamma13_luu;
    delete[] d1_gamma22_luu;    delete[] d2_gamma22_luu;    delete[] d3_gamma22_luu;
    delete[] d1_gamma23_luu;    delete[] d2_gamma23_luu;    delete[] d3_gamma23_luu;
    delete[] d1_gamma33_luu;    delete[] d2_gamma33_luu;    delete[] d3_gamma33_luu;
    delete[] K11_ll;
    delete[] K12_ll;
    delete[] K13_ll;
    delete[] K22_ll;
    delete[] K23_ll;
    delete[] K33_ll;
    // Tetrad Data:
    delete[] tetrad00_ul;    delete[] tetrad01_ul;    delete[] tetrad02_ul;    delete[] tetrad03_ul;
    delete[] tetrad10_ul;    delete[] tetrad11_ul;    delete[] tetrad12_ul;    delete[] tetrad13_ul;
    delete[] tetrad20_ul;    delete[] tetrad21_ul;    delete[] tetrad22_ul;    delete[] tetrad23_ul;
    delete[] tetrad30_ul;    delete[] tetrad31_ul;    delete[] tetrad32_ul;    delete[] tetrad33_ul;
}



std::string Metric::Name()
{
    exit_on_error("Metric virtual Method (Name) has been called!");
    return "";
}

// Initialization:
Tensor4x4 Metric::MetricFunction(const Coord& xyz)
{
    exit_on_error("Metric virtual Method (MetricFunction) has been called!");
    return Tensor4x4(0);
}

void Metric::InitializeMetricOnGrid()
{
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        Coord xyz = grid.xyz(i,j,k);
        Tensor4x4 g_ll = MetricFunction(xyz);
        Tensor4x4 g_uu = g_ll.Invert();
        int ijk = grid.Index(i,j,k);

        g00_ll[ijk] = g_ll[{0,0}];    g00_uu[ijk] = g_uu[{0,0}];
        g01_ll[ijk] = g_ll[{0,1}];    g01_uu[ijk] = g_uu[{0,1}];
        g02_ll[ijk] = g_ll[{0,2}];    g02_uu[ijk] = g_uu[{0,2}];
        g03_ll[ijk] = g_ll[{0,3}];    g03_uu[ijk] = g_uu[{0,3}];
        g11_ll[ijk] = g_ll[{1,1}];    g11_uu[ijk] = g_uu[{1,1}];
        g12_ll[ijk] = g_ll[{1,2}];    g12_uu[ijk] = g_uu[{1,2}];
        g13_ll[ijk] = g_ll[{1,3}];    g13_uu[ijk] = g_uu[{1,3}];
        g22_ll[ijk] = g_ll[{2,2}];    g22_uu[ijk] = g_uu[{2,2}]; 
        g23_ll[ijk] = g_ll[{2,3}];    g23_uu[ijk] = g_uu[{2,3}]; 
        g33_ll[ijk] = g_ll[{3,3}];    g33_uu[ijk] = g_uu[{3,3}]; 
    }
}

void Metric::InitializeBoostedTetradOnGrid()
{
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
        

        // Sub matrizes and inverses:
        Tensor3x3 A_ll = GetGamma_ll(ijk);
        Tensor3x3 A_uu = GetGamma_uu(ijk);
        Tensor2x2 B_ll(A_ll[{2,2}],A_ll[{2,3}], A_ll[{3,2}],A_ll[{3,3}]);
        Tensor2x2 B_uu = B_ll.Invert();

        // 'lapse' and 'shift' of sub matrizes:
        double alphaA = 1.0 / sqrt(A_uu[{1,1}]);
        double betaAx = alphaA * alphaA * A_uu[{1,2}];
        double betaAy = alphaA * alphaA * A_uu[{1,3}];
        
        double alphaB = 1.0 / sqrt(B_uu[{2,2}]);
        double betaB = alphaB * alphaB * B_uu[{2,3}];

        // Build eulerian observer like velocities:
        Tensor4 n4 = uEulObs(ijk);
        Tensor3 n3(1.0 / alphaA, betaAx / alphaA, betaAy / alphaA);
        Tensor2 n2(1.0 / alphaB, betaB  / alphaB);
        double  n1 =  1.0 / sqrt(B_ll[{3,3}]);

        // Construct tetrad:
        tetrad00_ul[ijk]=n4[0];    tetrad01_ul[ijk]=0;      tetrad02_ul[ijk] = 0;       tetrad03_ul[ijk] = 0;
        tetrad10_ul[ijk]=n4[1];    tetrad11_ul[ijk]=n3[1];  tetrad12_ul[ijk] = 0;       tetrad13_ul[ijk] = 0;
        tetrad20_ul[ijk]=n4[2];    tetrad21_ul[ijk]=n3[2];  tetrad22_ul[ijk] = n2[2];   tetrad23_ul[ijk] = 0;
        tetrad30_ul[ijk]=n4[3];    tetrad31_ul[ijk]=n3[3];  tetrad32_ul[ijk] = n2[3];   tetrad33_ul[ijk] = n1;
    }
}

template<int k>
Tensor4x4 Metric::MetricDeriv(const Coord& xyz)
{
    double dk = 1e-8;
    Tensor4x4 g_ll_m2 = MetricFunction(Coord{xyz[1]-2*dk*(k==1), xyz[2]-2*dk*(k==2), xyz[3]-2*dk*(k==3)});
    Tensor4x4 g_ll_m1 = MetricFunction(Coord{xyz[1]-1*dk*(k==1), xyz[2]-1*dk*(k==2), xyz[3]-1*dk*(k==3)});
    Tensor4x4 g_ll_p1 = MetricFunction(Coord{xyz[1]+1*dk*(k==1), xyz[2]+1*dk*(k==2), xyz[3]+1*dk*(k==3)});
    Tensor4x4 g_ll_p2 = MetricFunction(Coord{xyz[1]+2*dk*(k==1), xyz[2]+2*dk*(k==2), xyz[3]+2*dk*(k==3)});

    Tensor4x4 dk_g_ll;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
        dk_g_ll[{i,j}] = (1.0/12.0*g_ll_m2[{i,j}] - 1/12.0*g_ll_p2[{i,j}] - 2.0/3.0*g_ll_m1[{i,j}] + 2.0/3.0*g_ll_p1[{i,j}]) / dk;
    return dk_g_ll;
}

template<int k>
Tensor4x4 Metric::InverseMetricDeriv(const Coord& xyz)
{
    double dk = 1e-8;
    Tensor4x4 g_uu_m2 = MetricFunction(Coord{xyz[1]-2*dk*(k==1), xyz[2]-2*dk*(k==2), xyz[3]-2*dk*(k==3)}).Invert();
    Tensor4x4 g_uu_m1 = MetricFunction(Coord{xyz[1]-1*dk*(k==1), xyz[2]-1*dk*(k==2), xyz[3]-1*dk*(k==3)}).Invert();
    Tensor4x4 g_uu_p1 = MetricFunction(Coord{xyz[1]+1*dk*(k==1), xyz[2]+1*dk*(k==2), xyz[3]+1*dk*(k==3)}).Invert();
    Tensor4x4 g_uu_p2 = MetricFunction(Coord{xyz[1]+2*dk*(k==1), xyz[2]+2*dk*(k==2), xyz[3]+2*dk*(k==3)}).Invert();

    Tensor4x4 dk_g_uu;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
        dk_g_uu[{i,j}] = (1.0/12.0*g_uu_m2[{i,j}] - 1/12.0*g_uu_p2[{i,j}] - 2.0/3.0*g_uu_m1[{i,j}] + 2/3.0*g_uu_p1[{i,j}]) / dk;
    return dk_g_uu;
}

void Metric::InitializeMetricDerivativesOnGrid()
{
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        Coord xyz = grid.xyz(i,j,k);
        Tensor4x4 dt_gll(0.0);
        Tensor4x4 dx_gll = MetricDeriv<1>(xyz);
        Tensor4x4 dy_gll = MetricDeriv<2>(xyz);
        Tensor4x4 dz_gll = MetricDeriv<3>(xyz);
        Tensor4x4 dt_guu(0.0);
        Tensor4x4 dx_guu = InverseMetricDeriv<1>(xyz);
        Tensor4x4 dy_guu = InverseMetricDeriv<2>(xyz);
        Tensor4x4 dz_guu = InverseMetricDeriv<3>(xyz);

        int ijk = grid.Index(i,j,k);
        d0_g00_lll[ijk] = dt_gll[{0,0}];    d1_g00_lll[ijk] = dx_gll[{0,0}];     d2_g00_lll[ijk] = dy_gll[{0,0}];     d3_g00_lll[ijk] = dz_gll[{0,0}];
        d0_g01_lll[ijk] = dt_gll[{0,1}];    d1_g01_lll[ijk] = dx_gll[{0,1}];     d2_g01_lll[ijk] = dy_gll[{0,1}];     d3_g01_lll[ijk] = dz_gll[{0,1}];
        d0_g02_lll[ijk] = dt_gll[{0,2}];    d1_g02_lll[ijk] = dx_gll[{0,2}];     d2_g02_lll[ijk] = dy_gll[{0,2}];     d3_g02_lll[ijk] = dz_gll[{0,2}];
        d0_g03_lll[ijk] = dt_gll[{0,3}];    d1_g03_lll[ijk] = dx_gll[{0,3}];     d2_g03_lll[ijk] = dy_gll[{0,3}];     d3_g03_lll[ijk] = dz_gll[{0,3}];
        d0_g11_lll[ijk] = dt_gll[{1,1}];    d1_g11_lll[ijk] = dx_gll[{1,1}];     d2_g11_lll[ijk] = dy_gll[{1,1}];     d3_g11_lll[ijk] = dz_gll[{1,1}];
        d0_g12_lll[ijk] = dt_gll[{1,2}];    d1_g12_lll[ijk] = dx_gll[{1,2}];     d2_g12_lll[ijk] = dy_gll[{1,2}];     d3_g12_lll[ijk] = dz_gll[{1,2}];
        d0_g13_lll[ijk] = dt_gll[{1,3}];    d1_g13_lll[ijk] = dx_gll[{1,3}];     d2_g13_lll[ijk] = dy_gll[{1,3}];     d3_g13_lll[ijk] = dz_gll[{1,3}];
        d0_g22_lll[ijk] = dt_gll[{2,2}];    d1_g22_lll[ijk] = dx_gll[{2,2}];     d2_g22_lll[ijk] = dy_gll[{2,2}];     d3_g22_lll[ijk] = dz_gll[{2,2}];
        d0_g23_lll[ijk] = dt_gll[{2,3}];    d1_g23_lll[ijk] = dx_gll[{2,3}];     d2_g23_lll[ijk] = dy_gll[{2,3}];     d3_g23_lll[ijk] = dz_gll[{2,3}];
        d0_g33_lll[ijk] = dt_gll[{3,3}];    d1_g33_lll[ijk] = dx_gll[{3,3}];     d2_g33_lll[ijk] = dy_gll[{3,3}];     d3_g33_lll[ijk] = dz_gll[{3,3}];

        d0_g00_luu[ijk] = dt_guu[{0,0}];    d1_g00_luu[ijk] = dx_guu[{0,0}];     d2_g00_luu[ijk] = dy_guu[{0,0}];     d3_g00_luu[ijk] = dz_guu[{0,0}];
        d0_g01_luu[ijk] = dt_guu[{0,1}];    d1_g01_luu[ijk] = dx_guu[{0,1}];     d2_g01_luu[ijk] = dy_guu[{0,1}];     d3_g01_luu[ijk] = dz_guu[{0,1}];
        d0_g02_luu[ijk] = dt_guu[{0,2}];    d1_g02_luu[ijk] = dx_guu[{0,2}];     d2_g02_luu[ijk] = dy_guu[{0,2}];     d3_g02_luu[ijk] = dz_guu[{0,2}];
        d0_g03_luu[ijk] = dt_guu[{0,3}];    d1_g03_luu[ijk] = dx_guu[{0,3}];     d2_g03_luu[ijk] = dy_guu[{0,3}];     d3_g03_luu[ijk] = dz_guu[{0,3}];
        d0_g11_luu[ijk] = dt_guu[{1,1}];    d1_g11_luu[ijk] = dx_guu[{1,1}];     d2_g11_luu[ijk] = dy_guu[{1,1}];     d3_g11_luu[ijk] = dz_guu[{1,1}];
        d0_g12_luu[ijk] = dt_guu[{1,2}];    d1_g12_luu[ijk] = dx_guu[{1,2}];     d2_g12_luu[ijk] = dy_guu[{1,2}];     d3_g12_luu[ijk] = dz_guu[{1,2}];
        d0_g13_luu[ijk] = dt_guu[{1,3}];    d1_g13_luu[ijk] = dx_guu[{1,3}];     d2_g13_luu[ijk] = dy_guu[{1,3}];     d3_g13_luu[ijk] = dz_guu[{1,3}];
        d0_g22_luu[ijk] = dt_guu[{2,2}];    d1_g22_luu[ijk] = dx_guu[{2,2}];     d2_g22_luu[ijk] = dy_guu[{2,2}];     d3_g22_luu[ijk] = dz_guu[{2,2}];
        d0_g23_luu[ijk] = dt_guu[{2,3}];    d1_g23_luu[ijk] = dx_guu[{2,3}];     d2_g23_luu[ijk] = dy_guu[{2,3}];     d3_g23_luu[ijk] = dz_guu[{2,3}];
        d0_g33_luu[ijk] = dt_guu[{3,3}];    d1_g33_luu[ijk] = dx_guu[{3,3}];     d2_g33_luu[ijk] = dy_guu[{3,3}];     d3_g33_luu[ijk] = dz_guu[{3,3}];
    }
}

void Metric::InitializeAdmComponentsOnGrid()
{
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
        alpha[ijk] = 1.0/sqrt(-g00_uu[ijk]);
        beta1_u[ijk] = alpha[ijk] * alpha[ijk] * g01_uu[ijk];
        beta2_u[ijk] = alpha[ijk] * alpha[ijk] * g02_uu[ijk];
        beta3_u[ijk] = alpha[ijk] * alpha[ijk] * g03_uu[ijk];
        beta1_l[ijk] = g01_ll[ijk];
        beta2_l[ijk] = g02_ll[ijk];
        beta3_l[ijk] = g03_ll[ijk];
        gamma11_ll[ijk] = g11_ll[ijk];    gamma11_uu[ijk] = g11_uu[ijk] + beta1_u[ijk] * beta1_u[ijk] / (alpha[ijk] * alpha[ijk]);
        gamma12_ll[ijk] = g12_ll[ijk];    gamma12_uu[ijk] = g12_uu[ijk] + beta1_u[ijk] * beta2_u[ijk] / (alpha[ijk] * alpha[ijk]);
        gamma13_ll[ijk] = g13_ll[ijk];    gamma13_uu[ijk] = g13_uu[ijk] + beta1_u[ijk] * beta3_u[ijk] / (alpha[ijk] * alpha[ijk]);
        gamma22_ll[ijk] = g22_ll[ijk];    gamma22_uu[ijk] = g22_uu[ijk] + beta2_u[ijk] * beta2_u[ijk] / (alpha[ijk] * alpha[ijk]);
        gamma23_ll[ijk] = g23_ll[ijk];    gamma23_uu[ijk] = g23_uu[ijk] + beta2_u[ijk] * beta3_u[ijk] / (alpha[ijk] * alpha[ijk]);
        gamma33_ll[ijk] = g33_ll[ijk];    gamma33_uu[ijk] = g33_uu[ijk] + beta3_u[ijk] * beta3_u[ijk] / (alpha[ijk] * alpha[ijk]);

        d1_alpha_l[ijk] = d1_g00_luu[ijk] / (2.0 * pow(-g00_uu[ijk],3.0/2.0));
        d2_alpha_l[ijk] = d2_g00_luu[ijk] / (2.0 * pow(-g00_uu[ijk],3.0/2.0));
        d3_alpha_l[ijk] = d3_g00_luu[ijk] / (2.0 * pow(-g00_uu[ijk],3.0/2.0));
        d1_beta1_lu[ijk] = alpha[ijk] * alpha[ijk] * d1_g01_luu[ijk] + 2.0 * g01_uu[ijk] * alpha[ijk] * d1_alpha_l[ijk];
        d2_beta1_lu[ijk] = alpha[ijk] * alpha[ijk] * d2_g01_luu[ijk] + 2.0 * g01_uu[ijk] * alpha[ijk] * d2_alpha_l[ijk];
        d3_beta1_lu[ijk] = alpha[ijk] * alpha[ijk] * d3_g01_luu[ijk] + 2.0 * g01_uu[ijk] * alpha[ijk] * d3_alpha_l[ijk];
        d1_beta2_lu[ijk] = alpha[ijk] * alpha[ijk] * d1_g02_luu[ijk] + 2.0 * g02_uu[ijk] * alpha[ijk] * d1_alpha_l[ijk];
        d2_beta2_lu[ijk] = alpha[ijk] * alpha[ijk] * d2_g02_luu[ijk] + 2.0 * g02_uu[ijk] * alpha[ijk] * d2_alpha_l[ijk];
        d3_beta2_lu[ijk] = alpha[ijk] * alpha[ijk] * d3_g02_luu[ijk] + 2.0 * g02_uu[ijk] * alpha[ijk] * d3_alpha_l[ijk];
        d1_beta3_lu[ijk] = alpha[ijk] * alpha[ijk] * d1_g03_luu[ijk] + 2.0 * g03_uu[ijk] * alpha[ijk] * d1_alpha_l[ijk];
        d2_beta3_lu[ijk] = alpha[ijk] * alpha[ijk] * d2_g03_luu[ijk] + 2.0 * g03_uu[ijk] * alpha[ijk] * d2_alpha_l[ijk];
        d3_beta3_lu[ijk] = alpha[ijk] * alpha[ijk] * d3_g03_luu[ijk] + 2.0 * g03_uu[ijk] * alpha[ijk] * d3_alpha_l[ijk];
        d1_beta1_ll[ijk] = d1_g01_lll[ijk];    d1_beta2_ll[ijk] = d1_g02_lll[ijk];    d1_beta3_ll[ijk] = d1_g03_lll[ijk];
        d2_beta1_ll[ijk] = d2_g01_lll[ijk];    d2_beta2_ll[ijk] = d2_g02_lll[ijk];    d2_beta3_ll[ijk] = d2_g03_lll[ijk];
        d3_beta1_ll[ijk] = d3_g01_lll[ijk];    d3_beta2_ll[ijk] = d3_g02_lll[ijk];    d3_beta3_ll[ijk] = d3_g03_lll[ijk];
        d1_gamma11_lll[ijk] = d1_g11_lll[ijk];    d2_gamma11_lll[ijk] = d2_g11_lll[ijk];    d3_gamma11_lll[ijk] = d3_g11_lll[ijk];
        d1_gamma12_lll[ijk] = d1_g12_lll[ijk];    d2_gamma12_lll[ijk] = d2_g12_lll[ijk];    d3_gamma12_lll[ijk] = d3_g12_lll[ijk];
        d1_gamma13_lll[ijk] = d1_g13_lll[ijk];    d2_gamma13_lll[ijk] = d2_g13_lll[ijk];    d3_gamma13_lll[ijk] = d3_g13_lll[ijk];
        d1_gamma22_lll[ijk] = d1_g22_lll[ijk];    d2_gamma22_lll[ijk] = d2_g22_lll[ijk];    d3_gamma22_lll[ijk] = d3_g22_lll[ijk];
        d1_gamma23_lll[ijk] = d1_g23_lll[ijk];    d2_gamma23_lll[ijk] = d2_g23_lll[ijk];    d3_gamma23_lll[ijk] = d3_g23_lll[ijk];
        d1_gamma33_lll[ijk] = d1_g33_lll[ijk];    d2_gamma33_lll[ijk] = d2_g33_lll[ijk];    d3_gamma33_lll[ijk] = d3_g33_lll[ijk];
        d1_gamma11_luu[ijk] = d1_g11_luu[ijk] + (d1_beta1_lu[ijk] * beta1_u[ijk] + beta1_u[ijk] * d1_beta1_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta1_u[ijk] * d1_alpha_l[ijk]/pow(alpha[ijk],3);
        d1_gamma12_luu[ijk] = d1_g12_luu[ijk] + (d1_beta1_lu[ijk] * beta2_u[ijk] + beta1_u[ijk] * d1_beta2_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta2_u[ijk] * d1_alpha_l[ijk]/pow(alpha[ijk],3);
        d1_gamma13_luu[ijk] = d1_g13_luu[ijk] + (d1_beta1_lu[ijk] * beta3_u[ijk] + beta1_u[ijk] * d1_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta3_u[ijk] * d1_alpha_l[ijk]/pow(alpha[ijk],3);
        d1_gamma22_luu[ijk] = d1_g22_luu[ijk] + (d1_beta2_lu[ijk] * beta2_u[ijk] + beta2_u[ijk] * d1_beta2_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta2_u[ijk] * beta2_u[ijk] * d1_alpha_l[ijk]/pow(alpha[ijk],3);
        d1_gamma23_luu[ijk] = d1_g23_luu[ijk] + (d1_beta2_lu[ijk] * beta3_u[ijk] + beta2_u[ijk] * d1_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta2_u[ijk] * beta3_u[ijk] * d1_alpha_l[ijk]/pow(alpha[ijk],3);
        d1_gamma33_luu[ijk] = d1_g33_luu[ijk] + (d1_beta3_lu[ijk] * beta3_u[ijk] + beta3_u[ijk] * d1_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta3_u[ijk] * beta3_u[ijk] * d1_alpha_l[ijk]/pow(alpha[ijk],3);
        d2_gamma11_luu[ijk] = d2_g11_luu[ijk] + (d2_beta1_lu[ijk] * beta1_u[ijk] + beta1_u[ijk] * d2_beta1_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta1_u[ijk] * d2_alpha_l[ijk]/pow(alpha[ijk],3);
        d2_gamma12_luu[ijk] = d2_g12_luu[ijk] + (d2_beta1_lu[ijk] * beta2_u[ijk] + beta1_u[ijk] * d2_beta2_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta2_u[ijk] * d2_alpha_l[ijk]/pow(alpha[ijk],3);
        d2_gamma13_luu[ijk] = d2_g13_luu[ijk] + (d2_beta1_lu[ijk] * beta3_u[ijk] + beta1_u[ijk] * d2_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta3_u[ijk] * d2_alpha_l[ijk]/pow(alpha[ijk],3);
        d2_gamma22_luu[ijk] = d2_g22_luu[ijk] + (d2_beta2_lu[ijk] * beta2_u[ijk] + beta2_u[ijk] * d2_beta2_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta2_u[ijk] * beta2_u[ijk] * d2_alpha_l[ijk]/pow(alpha[ijk],3);
        d2_gamma23_luu[ijk] = d2_g23_luu[ijk] + (d2_beta2_lu[ijk] * beta3_u[ijk] + beta2_u[ijk] * d2_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta2_u[ijk] * beta3_u[ijk] * d2_alpha_l[ijk]/pow(alpha[ijk],3);
        d2_gamma33_luu[ijk] = d2_g33_luu[ijk] + (d2_beta3_lu[ijk] * beta3_u[ijk] + beta3_u[ijk] * d2_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta3_u[ijk] * beta3_u[ijk] * d2_alpha_l[ijk]/pow(alpha[ijk],3);  
        d3_gamma11_luu[ijk] = d3_g11_luu[ijk] + (d3_beta1_lu[ijk] * beta1_u[ijk] + beta1_u[ijk] * d3_beta1_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta1_u[ijk] * d3_alpha_l[ijk]/pow(alpha[ijk],3);
        d3_gamma12_luu[ijk] = d3_g12_luu[ijk] + (d3_beta1_lu[ijk] * beta2_u[ijk] + beta1_u[ijk] * d3_beta2_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta2_u[ijk] * d3_alpha_l[ijk]/pow(alpha[ijk],3);
        d3_gamma13_luu[ijk] = d3_g13_luu[ijk] + (d3_beta1_lu[ijk] * beta3_u[ijk] + beta1_u[ijk] * d3_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta1_u[ijk] * beta3_u[ijk] * d3_alpha_l[ijk]/pow(alpha[ijk],3);
        d3_gamma22_luu[ijk] = d3_g22_luu[ijk] + (d3_beta2_lu[ijk] * beta2_u[ijk] + beta2_u[ijk] * d3_beta2_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta2_u[ijk] * beta2_u[ijk] * d3_alpha_l[ijk]/pow(alpha[ijk],3);
        d3_gamma23_luu[ijk] = d3_g23_luu[ijk] + (d3_beta2_lu[ijk] * beta3_u[ijk] + beta2_u[ijk] * d3_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta2_u[ijk] * beta3_u[ijk] * d3_alpha_l[ijk]/pow(alpha[ijk],3);
        d3_gamma33_luu[ijk] = d3_g33_luu[ijk] + (d3_beta3_lu[ijk] * beta3_u[ijk] + beta3_u[ijk] * d3_beta3_lu[ijk]) / (alpha[ijk] * alpha[ijk]) - 2.0 * beta3_u[ijk] * beta3_u[ijk] * d3_alpha_l[ijk]/pow(alpha[ijk],3);

        Tensor3 beta_u = GetBeta_u(ijk);
        Tensor3x3 dl_beta_u = GetDerivBeta_lu(ijk);
        Tensor3x3 gamma_ll = GetGamma_ll(ijk);
        Tensor3x3 dt_gamma_ll(0.0);
        Tensor3x3x3 dl_gamma_ll = GetDerivGamma_lll(ijk);
        Tensor3x3 K_ll(0.0);
        for(int a=1; a<4; a++)
        for(int b=1; b<4; b++)
        {
            K_ll[{a,b}] = -dt_gamma_ll[{a,b}];
            for(int c=1; c<4; c++)
                K_ll[{a,b}] += 2.0 * gamma_ll[{a,c}] * dl_beta_u[{b,c}] + dl_gamma_ll[{c,a,b}] * beta_u[c];
            K_ll[{a,b}] /= 2.0 * alpha[ijk];
        }
        K11_ll[ijk] = K_ll[{1,1}];
        K12_ll[ijk] = K_ll[{1,2}];
        K13_ll[ijk] = K_ll[{1,3}];
        K22_ll[ijk] = K_ll[{2,2}];
        K23_ll[ijk] = K_ll[{2,3}];
        K33_ll[ijk] = K_ll[{3,3}];
    }
}

double Metric::InterpolateArrayTo_ijk(double* array, const Coord& ijk)
{
    int i0 = std::floor(ijk[1]);
    int j0 = std::floor(ijk[2]);
    int k0 = std::floor(ijk[3]);
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;

    return TrilinearInterpolation
    (ijk[1] - i0, ijk[2] - j0, ijk[3] - k0,
    array[grid.Index(i0,j0,k0)], array[grid.Index(i0,j0,k1)], array[grid.Index(i0,j1,k0)], array[grid.Index(i0,j1,k1)],
    array[grid.Index(i1,j0,k0)], array[grid.Index(i1,j0,k1)], array[grid.Index(i1,j1,k0)], array[grid.Index(i1,j1,k1)]);
}
double Metric::InterpolateArrayTo_ijk(double* array, double i, double j, double k)
{
    int i0 = std::floor(i);
    int j0 = std::floor(j);
    int k0 = std::floor(k);
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;

    return TrilinearInterpolation
    (i - i0, j - j0, k - k0,
    array[grid.Index(i0,j0,k0)], array[grid.Index(i0,j0,k1)], array[grid.Index(i0,j1,k0)], array[grid.Index(i0,j1,k1)],
    array[grid.Index(i1,j0,k0)], array[grid.Index(i1,j0,k1)], array[grid.Index(i1,j1,k0)], array[grid.Index(i1,j1,k1)]);
}

bool Metric::InsideBH(const Coord& xyz)
{
    exit_on_error("Metric virtual Method (InsideBH(xyz)) has been called!");
    return false;
}



// Tensor getters:
Tensor4 Metric::uEulObs(int ijk)
{
    double alpha = GetAlpha(ijk);
    Tensor3 beta_u = GetBeta_u(ijk);
    return Tensor4(1.0/alpha, -beta_u[1]/alpha, -beta_u[2]/alpha, -beta_u[3]/alpha);
}
Tensor4 Metric::uEulObs(const Coord& xyz)
{
    double alpha = GetAlpha(xyz);
    Tensor3 beta_u = GetBeta_u(xyz);
    if(grid.OutsideDomain(xyz))
        return Tensor4(1,0,0,0);
    else
        return Tensor4 (1.0/alpha, -beta_u[1]/alpha, -beta_u[2]/alpha, -beta_u[3]/alpha);
}
Tensor4x4 Metric::GetMetric_ll(int ijk)
{
    return Tensor4x4
    (g00_ll[ijk], g01_ll[ijk], g02_ll[ijk], g03_ll[ijk],
     g01_ll[ijk], g11_ll[ijk], g12_ll[ijk], g13_ll[ijk],
     g02_ll[ijk], g12_ll[ijk], g22_ll[ijk], g23_ll[ijk],
     g03_ll[ijk], g13_ll[ijk], g23_ll[ijk], g33_ll[ijk]);
}
Tensor4x4 Metric::GetMetric_ll(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor4x4(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor4x4 g_ll;
        g_ll[{0,0}] = InterpolateArrayTo_ijk(g00_ll,ijk);
        g_ll[{1,1}] = InterpolateArrayTo_ijk(g11_ll,ijk);
        g_ll[{2,2}] = InterpolateArrayTo_ijk(g22_ll,ijk);
        g_ll[{3,3}] = InterpolateArrayTo_ijk(g33_ll,ijk);
        g_ll[{0,1}] = g_ll[{1,0}] = InterpolateArrayTo_ijk(g01_ll,ijk);
        g_ll[{0,2}] = g_ll[{2,0}] = InterpolateArrayTo_ijk(g02_ll,ijk);
        g_ll[{0,3}] = g_ll[{3,0}] = InterpolateArrayTo_ijk(g03_ll,ijk);
        g_ll[{1,2}] = g_ll[{2,1}] = InterpolateArrayTo_ijk(g12_ll,ijk);
        g_ll[{1,3}] = g_ll[{3,1}] = InterpolateArrayTo_ijk(g13_ll,ijk);
        g_ll[{2,3}] = g_ll[{3,2}] = InterpolateArrayTo_ijk(g23_ll,ijk);
        return g_ll;
    }
}
Tensor4x4 Metric::GetMetric_uu(int ijk)
{
    return Tensor4x4
    (g00_uu[ijk], g01_uu[ijk], g02_uu[ijk], g03_uu[ijk],
     g01_uu[ijk], g11_uu[ijk], g12_uu[ijk], g13_uu[ijk],
     g02_uu[ijk], g12_uu[ijk], g22_uu[ijk], g23_uu[ijk],
     g03_uu[ijk], g13_uu[ijk], g23_uu[ijk], g33_uu[ijk]);
}
Tensor4x4 Metric::GetMetric_uu(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor4x4(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor4x4 g_uu;
        g_uu[{0,0}] = InterpolateArrayTo_ijk(g00_uu,ijk);
        g_uu[{1,1}] = InterpolateArrayTo_ijk(g11_uu,ijk);
        g_uu[{2,2}] = InterpolateArrayTo_ijk(g22_uu,ijk);
        g_uu[{3,3}] = InterpolateArrayTo_ijk(g33_uu,ijk);
        g_uu[{0,1}] = g_uu[{1,0}] = InterpolateArrayTo_ijk(g01_uu,ijk);
        g_uu[{0,2}] = g_uu[{2,0}] = InterpolateArrayTo_ijk(g02_uu,ijk);
        g_uu[{0,3}] = g_uu[{3,0}] = InterpolateArrayTo_ijk(g03_uu,ijk);
        g_uu[{1,2}] = g_uu[{2,1}] = InterpolateArrayTo_ijk(g12_uu,ijk);
        g_uu[{1,3}] = g_uu[{3,1}] = InterpolateArrayTo_ijk(g13_uu,ijk);
        g_uu[{2,3}] = g_uu[{3,2}] = InterpolateArrayTo_ijk(g23_uu,ijk);
        return g_uu;
    }
}

Tensor4x4 Metric::GetMinkowskiMetric_ll(int ijk)
{
    return Tensor4x4(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
}
Tensor4x4 Metric::GetMinkowskiMetric_ll(const Coord& xyz)
{
    return Tensor4x4(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
}
Tensor4x4 Metric::GetMinkowskiMetric_uu(int ijk)
{
    return Tensor4x4(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
}
Tensor4x4 Metric::GetMinkowskiMetric_uu(const Coord& xyz)
{
    return Tensor4x4(-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
}

Tensor4x4x4 Metric::GetDerivMetric_lll(int ijk)
{
    return Tensor4x4x4
    (d0_g00_lll[ijk], d0_g01_lll[ijk], d0_g02_lll[ijk], d0_g03_lll[ijk],
     d0_g01_lll[ijk], d0_g11_lll[ijk], d0_g12_lll[ijk], d0_g13_lll[ijk],
     d0_g02_lll[ijk], d0_g12_lll[ijk], d0_g22_lll[ijk], d0_g23_lll[ijk],
     d0_g03_lll[ijk], d0_g13_lll[ijk], d0_g23_lll[ijk], d0_g33_lll[ijk],

     d1_g00_lll[ijk], d1_g01_lll[ijk], d1_g02_lll[ijk], d1_g03_lll[ijk],
     d1_g01_lll[ijk], d1_g11_lll[ijk], d1_g12_lll[ijk], d1_g13_lll[ijk],
     d1_g02_lll[ijk], d1_g12_lll[ijk], d1_g22_lll[ijk], d1_g23_lll[ijk],
     d1_g03_lll[ijk], d1_g13_lll[ijk], d1_g23_lll[ijk], d1_g33_lll[ijk],

     d2_g00_lll[ijk], d2_g01_lll[ijk], d2_g02_lll[ijk], d2_g03_lll[ijk],
     d2_g01_lll[ijk], d2_g11_lll[ijk], d2_g12_lll[ijk], d2_g13_lll[ijk],
     d2_g02_lll[ijk], d2_g12_lll[ijk], d2_g22_lll[ijk], d2_g23_lll[ijk],
     d2_g03_lll[ijk], d2_g13_lll[ijk], d2_g23_lll[ijk], d2_g33_lll[ijk],

     d3_g00_lll[ijk], d3_g01_lll[ijk], d3_g02_lll[ijk], d3_g03_lll[ijk],
     d3_g01_lll[ijk], d3_g11_lll[ijk], d3_g12_lll[ijk], d3_g13_lll[ijk],
     d3_g02_lll[ijk], d3_g12_lll[ijk], d3_g22_lll[ijk], d3_g23_lll[ijk],
     d3_g03_lll[ijk], d3_g13_lll[ijk], d3_g23_lll[ijk], d3_g33_lll[ijk]);
}

Tensor4x4x4 Metric::GetDerivMetric_lll(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor4x4x4(0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor4x4x4 dl_g_ll;
        dl_g_ll[{0,0,0}] = InterpolateArrayTo_ijk(d0_g00_lll,ijk);
        dl_g_ll[{0,1,1}] = InterpolateArrayTo_ijk(d0_g11_lll,ijk);
        dl_g_ll[{0,2,2}] = InterpolateArrayTo_ijk(d0_g22_lll,ijk);
        dl_g_ll[{0,3,3}] = InterpolateArrayTo_ijk(d0_g33_lll,ijk);
        dl_g_ll[{0,0,1}] = dl_g_ll[{0,1,0}] = InterpolateArrayTo_ijk(d0_g01_lll,ijk);
        dl_g_ll[{0,0,2}] = dl_g_ll[{0,2,0}] = InterpolateArrayTo_ijk(d0_g02_lll,ijk);
        dl_g_ll[{0,0,3}] = dl_g_ll[{0,3,0}] = InterpolateArrayTo_ijk(d0_g03_lll,ijk);
        dl_g_ll[{0,1,2}] = dl_g_ll[{0,2,1}] = InterpolateArrayTo_ijk(d0_g12_lll,ijk);
        dl_g_ll[{0,1,3}] = dl_g_ll[{0,3,1}] = InterpolateArrayTo_ijk(d0_g13_lll,ijk);
        dl_g_ll[{0,2,3}] = dl_g_ll[{0,3,2}] = InterpolateArrayTo_ijk(d0_g23_lll,ijk);

        dl_g_ll[{1,0,0}] = InterpolateArrayTo_ijk(d1_g00_lll,ijk);
        dl_g_ll[{1,1,1}] = InterpolateArrayTo_ijk(d1_g11_lll,ijk);
        dl_g_ll[{1,2,2}] = InterpolateArrayTo_ijk(d1_g22_lll,ijk);
        dl_g_ll[{1,3,3}] = InterpolateArrayTo_ijk(d1_g33_lll,ijk);
        dl_g_ll[{1,0,1}] = dl_g_ll[{1,1,0}] = InterpolateArrayTo_ijk(d1_g01_lll,ijk);
        dl_g_ll[{1,0,2}] = dl_g_ll[{1,2,0}] = InterpolateArrayTo_ijk(d1_g02_lll,ijk);
        dl_g_ll[{1,0,3}] = dl_g_ll[{1,3,0}] = InterpolateArrayTo_ijk(d1_g03_lll,ijk);
        dl_g_ll[{1,1,2}] = dl_g_ll[{1,2,1}] = InterpolateArrayTo_ijk(d1_g12_lll,ijk);
        dl_g_ll[{1,1,3}] = dl_g_ll[{1,3,1}] = InterpolateArrayTo_ijk(d1_g13_lll,ijk);
        dl_g_ll[{1,2,3}] = dl_g_ll[{1,3,2}] = InterpolateArrayTo_ijk(d1_g23_lll,ijk);

        dl_g_ll[{2,0,0}] = InterpolateArrayTo_ijk(d2_g00_lll,ijk);
        dl_g_ll[{2,1,1}] = InterpolateArrayTo_ijk(d2_g11_lll,ijk);
        dl_g_ll[{2,2,2}] = InterpolateArrayTo_ijk(d2_g22_lll,ijk);
        dl_g_ll[{2,3,3}] = InterpolateArrayTo_ijk(d2_g33_lll,ijk);
        dl_g_ll[{2,0,1}] = dl_g_ll[{2,1,0}] = InterpolateArrayTo_ijk(d2_g01_lll,ijk);
        dl_g_ll[{2,0,2}] = dl_g_ll[{2,2,0}] = InterpolateArrayTo_ijk(d2_g02_lll,ijk);
        dl_g_ll[{2,0,3}] = dl_g_ll[{2,3,0}] = InterpolateArrayTo_ijk(d2_g03_lll,ijk);
        dl_g_ll[{2,1,2}] = dl_g_ll[{2,2,1}] = InterpolateArrayTo_ijk(d2_g12_lll,ijk);
        dl_g_ll[{2,1,3}] = dl_g_ll[{2,3,1}] = InterpolateArrayTo_ijk(d2_g13_lll,ijk);
        dl_g_ll[{2,2,3}] = dl_g_ll[{2,3,2}] = InterpolateArrayTo_ijk(d2_g23_lll,ijk);

        dl_g_ll[{3,0,0}] = InterpolateArrayTo_ijk(d3_g00_lll,ijk);
        dl_g_ll[{3,1,1}] = InterpolateArrayTo_ijk(d3_g11_lll,ijk);
        dl_g_ll[{3,2,2}] = InterpolateArrayTo_ijk(d3_g22_lll,ijk);
        dl_g_ll[{3,3,3}] = InterpolateArrayTo_ijk(d3_g33_lll,ijk);
        dl_g_ll[{3,0,1}] = dl_g_ll[{3,1,0}] = InterpolateArrayTo_ijk(d3_g01_lll,ijk);
        dl_g_ll[{3,0,2}] = dl_g_ll[{3,2,0}] = InterpolateArrayTo_ijk(d3_g02_lll,ijk);
        dl_g_ll[{3,0,3}] = dl_g_ll[{3,3,0}] = InterpolateArrayTo_ijk(d3_g03_lll,ijk);
        dl_g_ll[{3,1,2}] = dl_g_ll[{3,2,1}] = InterpolateArrayTo_ijk(d3_g12_lll,ijk);
        dl_g_ll[{3,1,3}] = dl_g_ll[{3,3,1}] = InterpolateArrayTo_ijk(d3_g13_lll,ijk);
        dl_g_ll[{3,2,3}] = dl_g_ll[{3,3,2}] = InterpolateArrayTo_ijk(d3_g23_lll,ijk);
        return dl_g_ll;
    }
}



Tensor4x4x4 Metric::GetDerivMetric_luu(int ijk)
{
    return Tensor4x4x4
    (d0_g00_luu[ijk], d0_g01_luu[ijk], d0_g02_luu[ijk], d0_g03_luu[ijk],
     d0_g01_luu[ijk], d0_g11_luu[ijk], d0_g12_luu[ijk], d0_g13_luu[ijk],
     d0_g02_luu[ijk], d0_g12_luu[ijk], d0_g22_luu[ijk], d0_g23_luu[ijk],
     d0_g03_luu[ijk], d0_g13_luu[ijk], d0_g23_luu[ijk], d0_g33_luu[ijk],

     d1_g00_luu[ijk], d1_g01_luu[ijk], d1_g02_luu[ijk], d1_g03_luu[ijk],
     d1_g01_luu[ijk], d1_g11_luu[ijk], d1_g12_luu[ijk], d1_g13_luu[ijk],
     d1_g02_luu[ijk], d1_g12_luu[ijk], d1_g22_luu[ijk], d1_g23_luu[ijk],
     d1_g03_luu[ijk], d1_g13_luu[ijk], d1_g23_luu[ijk], d1_g33_luu[ijk],

     d2_g00_luu[ijk], d2_g01_luu[ijk], d2_g02_luu[ijk], d2_g03_luu[ijk],
     d2_g01_luu[ijk], d2_g11_luu[ijk], d2_g12_luu[ijk], d2_g13_luu[ijk],
     d2_g02_luu[ijk], d2_g12_luu[ijk], d2_g22_luu[ijk], d2_g23_luu[ijk],
     d2_g03_luu[ijk], d2_g13_luu[ijk], d2_g23_luu[ijk], d2_g33_luu[ijk],

     d3_g00_luu[ijk], d3_g01_luu[ijk], d3_g02_luu[ijk], d3_g03_luu[ijk],
     d3_g01_luu[ijk], d3_g11_luu[ijk], d3_g12_luu[ijk], d3_g13_luu[ijk],
     d3_g02_luu[ijk], d3_g12_luu[ijk], d3_g22_luu[ijk], d3_g23_luu[ijk],
     d3_g03_luu[ijk], d3_g13_luu[ijk], d3_g23_luu[ijk], d3_g33_luu[ijk]);
}

Tensor4x4x4 Metric::GetDerivMetric_luu(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor4x4x4(0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor4x4x4 dl_g_uu;
        dl_g_uu[{0,0,0}] = InterpolateArrayTo_ijk(d0_g00_luu,ijk);
        dl_g_uu[{0,1,1}] = InterpolateArrayTo_ijk(d0_g11_luu,ijk);
        dl_g_uu[{0,2,2}] = InterpolateArrayTo_ijk(d0_g22_luu,ijk);
        dl_g_uu[{0,3,3}] = InterpolateArrayTo_ijk(d0_g33_luu,ijk);
        dl_g_uu[{0,0,1}] = dl_g_uu[{0,1,0}] = InterpolateArrayTo_ijk(d0_g01_luu,ijk);
        dl_g_uu[{0,0,2}] = dl_g_uu[{0,2,0}] = InterpolateArrayTo_ijk(d0_g02_luu,ijk);
        dl_g_uu[{0,0,3}] = dl_g_uu[{0,3,0}] = InterpolateArrayTo_ijk(d0_g03_luu,ijk);
        dl_g_uu[{0,1,2}] = dl_g_uu[{0,2,1}] = InterpolateArrayTo_ijk(d0_g12_luu,ijk);
        dl_g_uu[{0,1,3}] = dl_g_uu[{0,3,1}] = InterpolateArrayTo_ijk(d0_g13_luu,ijk);
        dl_g_uu[{0,2,3}] = dl_g_uu[{0,3,2}] = InterpolateArrayTo_ijk(d0_g23_luu,ijk);

        dl_g_uu[{1,0,0}] = InterpolateArrayTo_ijk(d1_g00_luu,ijk);
        dl_g_uu[{1,1,1}] = InterpolateArrayTo_ijk(d1_g11_luu,ijk);
        dl_g_uu[{1,2,2}] = InterpolateArrayTo_ijk(d1_g22_luu,ijk);
        dl_g_uu[{1,3,3}] = InterpolateArrayTo_ijk(d1_g33_luu,ijk);
        dl_g_uu[{1,0,1}] = dl_g_uu[{1,1,0}] = InterpolateArrayTo_ijk(d1_g01_luu,ijk);
        dl_g_uu[{1,0,2}] = dl_g_uu[{1,2,0}] = InterpolateArrayTo_ijk(d1_g02_luu,ijk);
        dl_g_uu[{1,0,3}] = dl_g_uu[{1,3,0}] = InterpolateArrayTo_ijk(d1_g03_luu,ijk);
        dl_g_uu[{1,1,2}] = dl_g_uu[{1,2,1}] = InterpolateArrayTo_ijk(d1_g12_luu,ijk);
        dl_g_uu[{1,1,3}] = dl_g_uu[{1,3,1}] = InterpolateArrayTo_ijk(d1_g13_luu,ijk);
        dl_g_uu[{1,2,3}] = dl_g_uu[{1,3,2}] = InterpolateArrayTo_ijk(d1_g23_luu,ijk);

        dl_g_uu[{2,0,0}] = InterpolateArrayTo_ijk(d2_g00_luu,ijk);
        dl_g_uu[{2,1,1}] = InterpolateArrayTo_ijk(d2_g11_luu,ijk);
        dl_g_uu[{2,2,2}] = InterpolateArrayTo_ijk(d2_g22_luu,ijk);
        dl_g_uu[{2,3,3}] = InterpolateArrayTo_ijk(d2_g33_luu,ijk);
        dl_g_uu[{2,0,1}] = dl_g_uu[{2,1,0}] = InterpolateArrayTo_ijk(d2_g01_luu,ijk);
        dl_g_uu[{2,0,2}] = dl_g_uu[{2,2,0}] = InterpolateArrayTo_ijk(d2_g02_luu,ijk);
        dl_g_uu[{2,0,3}] = dl_g_uu[{2,3,0}] = InterpolateArrayTo_ijk(d2_g03_luu,ijk);
        dl_g_uu[{2,1,2}] = dl_g_uu[{2,2,1}] = InterpolateArrayTo_ijk(d2_g12_luu,ijk);
        dl_g_uu[{2,1,3}] = dl_g_uu[{2,3,1}] = InterpolateArrayTo_ijk(d2_g13_luu,ijk);
        dl_g_uu[{2,2,3}] = dl_g_uu[{2,3,2}] = InterpolateArrayTo_ijk(d2_g23_luu,ijk);

        dl_g_uu[{3,0,0}] = InterpolateArrayTo_ijk(d3_g00_luu,ijk);
        dl_g_uu[{3,1,1}] = InterpolateArrayTo_ijk(d3_g11_luu,ijk);
        dl_g_uu[{3,2,2}] = InterpolateArrayTo_ijk(d3_g22_luu,ijk);
        dl_g_uu[{3,3,3}] = InterpolateArrayTo_ijk(d3_g33_luu,ijk);
        dl_g_uu[{3,0,1}] = dl_g_uu[{3,1,0}] = InterpolateArrayTo_ijk(d3_g01_luu,ijk);
        dl_g_uu[{3,0,2}] = dl_g_uu[{3,2,0}] = InterpolateArrayTo_ijk(d3_g02_luu,ijk);
        dl_g_uu[{3,0,3}] = dl_g_uu[{3,3,0}] = InterpolateArrayTo_ijk(d3_g03_luu,ijk);
        dl_g_uu[{3,1,2}] = dl_g_uu[{3,2,1}] = InterpolateArrayTo_ijk(d3_g12_luu,ijk);
        dl_g_uu[{3,1,3}] = dl_g_uu[{3,3,1}] = InterpolateArrayTo_ijk(d3_g13_luu,ijk);
        dl_g_uu[{3,2,3}] = dl_g_uu[{3,3,2}] = InterpolateArrayTo_ijk(d3_g23_luu,ijk);
        return dl_g_uu;
    }
}


Tensor4x4 Metric::GetTetrad(int ijk)
{
    return Tensor4x4
    (tetrad00_ul[ijk], tetrad01_ul[ijk], tetrad02_ul[ijk], tetrad03_ul[ijk],
     tetrad10_ul[ijk], tetrad11_ul[ijk], tetrad12_ul[ijk], tetrad13_ul[ijk],
     tetrad20_ul[ijk], tetrad21_ul[ijk], tetrad22_ul[ijk], tetrad23_ul[ijk],
     tetrad30_ul[ijk], tetrad31_ul[ijk], tetrad32_ul[ijk], tetrad33_ul[ijk]);
}

Tensor4x4 Metric::GetTetrad(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor4x4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
    else
    {
        Coord ijk = grid.ijk(xyz);
        return Tensor4x4
        (InterpolateArrayTo_ijk(tetrad00_ul,ijk), InterpolateArrayTo_ijk(tetrad01_ul,ijk), InterpolateArrayTo_ijk(tetrad02_ul,ijk), InterpolateArrayTo_ijk(tetrad03_ul,ijk),
         InterpolateArrayTo_ijk(tetrad10_ul,ijk), InterpolateArrayTo_ijk(tetrad11_ul,ijk), InterpolateArrayTo_ijk(tetrad12_ul,ijk), InterpolateArrayTo_ijk(tetrad13_ul,ijk),
         InterpolateArrayTo_ijk(tetrad20_ul,ijk), InterpolateArrayTo_ijk(tetrad21_ul,ijk), InterpolateArrayTo_ijk(tetrad22_ul,ijk), InterpolateArrayTo_ijk(tetrad23_ul,ijk),
         InterpolateArrayTo_ijk(tetrad30_ul,ijk), InterpolateArrayTo_ijk(tetrad31_ul,ijk), InterpolateArrayTo_ijk(tetrad32_ul,ijk), InterpolateArrayTo_ijk(tetrad33_ul,ijk));
    }
}

// ADM getters:
double Metric::GetAlpha(int ijk)
{
    return alpha[ijk];
}
double Metric::GetAlpha(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return 1.0;
    else
    {
        Coord ijk = grid.ijk(xyz);
        return InterpolateArrayTo_ijk(alpha,ijk);
    }
}

Tensor3 Metric::GetBeta_u(int ijk)
{
    return Tensor3(beta1_u[ijk],beta2_u[ijk],beta3_u[ijk]);
}
Tensor3 Metric::GetBeta_u(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        return Tensor3(InterpolateArrayTo_ijk(beta1_u,ijk), InterpolateArrayTo_ijk(beta2_u,ijk), InterpolateArrayTo_ijk(beta3_u,ijk));
    }
}

Tensor3 Metric::GetBeta_l(int ijk)
{
    return Tensor3(beta1_l[ijk],beta2_l[ijk],beta3_l[ijk]);
}
Tensor3 Metric::GetBeta_l(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        return Tensor3(InterpolateArrayTo_ijk(beta1_l,ijk), InterpolateArrayTo_ijk(beta2_l,ijk), InterpolateArrayTo_ijk(beta3_l,ijk));
    }
}

Tensor3x3 Metric::GetGamma_ll(int ijk)
{
    return Tensor3x3
    (gamma11_ll[ijk],gamma12_ll[ijk],gamma13_ll[ijk],
     gamma12_ll[ijk],gamma22_ll[ijk],gamma23_ll[ijk],
     gamma13_ll[ijk],gamma23_ll[ijk],gamma33_ll[ijk]);
}
Tensor3x3 Metric::GetGamma_ll(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3(1,0,0, 0,1,0, 0,0,1);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor3x3 gamma_ll;
        gamma_ll[{1,1}] = InterpolateArrayTo_ijk(gamma11_ll,ijk);
        gamma_ll[{2,2}] = InterpolateArrayTo_ijk(gamma22_ll,ijk);
        gamma_ll[{3,3}] = InterpolateArrayTo_ijk(gamma33_ll,ijk);
        gamma_ll[{1,2}] = gamma_ll[{2,1}] = InterpolateArrayTo_ijk(gamma12_ll,ijk);
        gamma_ll[{1,3}] = gamma_ll[{3,1}] = InterpolateArrayTo_ijk(gamma13_ll,ijk);
        gamma_ll[{2,3}] = gamma_ll[{3,2}] = InterpolateArrayTo_ijk(gamma23_ll,ijk);
        return gamma_ll;
    }
}

Tensor3x3 Metric::GetGamma_uu(int ijk)
{
    return Tensor3x3
    (gamma11_uu[ijk],gamma12_uu[ijk],gamma13_uu[ijk],
     gamma12_uu[ijk],gamma22_uu[ijk],gamma23_uu[ijk],
     gamma13_uu[ijk],gamma23_uu[ijk],gamma33_uu[ijk]);
}
Tensor3x3 Metric::GetGamma_uu(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3(1,0,0, 0,1,0, 0,0,1);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor3x3 gamma_uu;
        gamma_uu[{1,1}] = InterpolateArrayTo_ijk(gamma11_uu,ijk);
        gamma_uu[{2,2}] = InterpolateArrayTo_ijk(gamma22_uu,ijk);
        gamma_uu[{3,3}] = InterpolateArrayTo_ijk(gamma33_uu,ijk);
        gamma_uu[{1,2}] = gamma_uu[{2,1}] = InterpolateArrayTo_ijk(gamma12_uu,ijk);
        gamma_uu[{1,3}] = gamma_uu[{3,1}] = InterpolateArrayTo_ijk(gamma13_uu,ijk);
        gamma_uu[{2,3}] = gamma_uu[{3,2}] = InterpolateArrayTo_ijk(gamma23_uu,ijk);
        return gamma_uu;
    }
}

Tensor3x3 Metric::GetMinkowskiGamma_ll(int ijk)
{
    return Tensor3x3(1,0,0, 0,1,0, 0,0,1);
}
Tensor3x3 Metric::GetMinkowskiGamma_ll(const Coord& xyz)
{
    return Tensor3x3(1,0,0, 0,1,0, 0,0,1);
}

Tensor3x3 Metric::GetMinkowskiGamma_uu(int ijk)
{
    return Tensor3x3(1,0,0, 0,1,0, 0,0,1);
}
Tensor3x3 Metric::GetMinkowskiGamma_uu(const Coord& xyz)
{
    return Tensor3x3(1,0,0, 0,1,0, 0,0,1);
}

Tensor3 Metric::GetDerivAlpha_l(int ijk)
{
    return Tensor3(d1_alpha_l[ijk],d2_alpha_l[ijk],d3_alpha_l[ijk]);
}
Tensor3 Metric::GetDerivAlpha_l(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        return Tensor3(InterpolateArrayTo_ijk(d1_alpha_l,ijk), InterpolateArrayTo_ijk(d2_alpha_l,ijk), InterpolateArrayTo_ijk(d3_alpha_l,ijk));
    }
}

Tensor3x3 Metric::GetDerivBeta_lu(int ijk)
{
    return Tensor3x3
    (d1_beta1_lu[ijk],d1_beta2_lu[ijk],d1_beta3_lu[ijk],
     d2_beta1_lu[ijk],d2_beta2_lu[ijk],d2_beta3_lu[ijk],
     d3_beta1_lu[ijk],d3_beta2_lu[ijk],d3_beta3_lu[ijk]);
}
Tensor3x3 Metric::GetDerivBeta_lu(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        return Tensor3x3
        (InterpolateArrayTo_ijk(d1_beta1_lu,ijk), InterpolateArrayTo_ijk(d1_beta2_lu,ijk), InterpolateArrayTo_ijk(d1_beta3_lu,ijk),
         InterpolateArrayTo_ijk(d2_beta1_lu,ijk), InterpolateArrayTo_ijk(d2_beta2_lu,ijk), InterpolateArrayTo_ijk(d2_beta3_lu,ijk),
         InterpolateArrayTo_ijk(d3_beta1_lu,ijk), InterpolateArrayTo_ijk(d3_beta2_lu,ijk), InterpolateArrayTo_ijk(d3_beta3_lu,ijk));
    }
}

Tensor3x3 Metric::GetDerivBeta_ll(int ijk)
{
    return Tensor3x3
    (d1_beta1_ll[ijk],d1_beta2_ll[ijk],d1_beta3_ll[ijk],
     d2_beta1_ll[ijk],d2_beta2_ll[ijk],d2_beta3_ll[ijk],
     d3_beta1_ll[ijk],d3_beta2_ll[ijk],d3_beta3_ll[ijk]);
}
Tensor3x3 Metric::GetDerivBeta_ll(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        return Tensor3x3
        (InterpolateArrayTo_ijk(d1_beta1_ll,ijk), InterpolateArrayTo_ijk(d1_beta2_ll,ijk), InterpolateArrayTo_ijk(d1_beta3_ll,ijk),
         InterpolateArrayTo_ijk(d2_beta1_ll,ijk), InterpolateArrayTo_ijk(d2_beta2_ll,ijk), InterpolateArrayTo_ijk(d2_beta3_ll,ijk),
         InterpolateArrayTo_ijk(d3_beta1_ll,ijk), InterpolateArrayTo_ijk(d3_beta2_ll,ijk), InterpolateArrayTo_ijk(d3_beta3_ll,ijk));
    }
}

Tensor3x3x3 Metric::GetDerivGamma_lll(int ijk)
{
    return Tensor3x3x3
    (d1_gamma11_lll[ijk], d1_gamma12_lll[ijk], d1_gamma13_lll[ijk],
     d1_gamma12_lll[ijk], d1_gamma22_lll[ijk], d1_gamma23_lll[ijk],
     d1_gamma13_lll[ijk], d1_gamma23_lll[ijk], d1_gamma33_lll[ijk],
     d2_gamma11_lll[ijk], d2_gamma12_lll[ijk], d2_gamma13_lll[ijk],
     d2_gamma12_lll[ijk], d2_gamma22_lll[ijk], d2_gamma23_lll[ijk],
     d2_gamma13_lll[ijk], d2_gamma23_lll[ijk], d2_gamma33_lll[ijk],
     d3_gamma11_lll[ijk], d3_gamma12_lll[ijk], d3_gamma13_lll[ijk],
     d3_gamma12_lll[ijk], d3_gamma22_lll[ijk], d3_gamma23_lll[ijk],
     d3_gamma13_lll[ijk], d3_gamma23_lll[ijk], d3_gamma33_lll[ijk]);
}
Tensor3x3x3 Metric::GetDerivGamma_lll(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3x3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor3x3x3 dGamma_lll;
        dGamma_lll[{1,1,1}] = InterpolateArrayTo_ijk(d1_gamma11_lll,ijk);
        dGamma_lll[{1,2,2}] = InterpolateArrayTo_ijk(d1_gamma22_lll,ijk);
        dGamma_lll[{1,3,3}] = InterpolateArrayTo_ijk(d1_gamma33_lll,ijk);
        dGamma_lll[{1,1,2}] = dGamma_lll[{1,2,1}] = InterpolateArrayTo_ijk(d1_gamma12_lll,ijk);
        dGamma_lll[{1,1,3}] = dGamma_lll[{1,3,1}] = InterpolateArrayTo_ijk(d1_gamma13_lll,ijk);
        dGamma_lll[{1,2,3}] = dGamma_lll[{1,3,2}] = InterpolateArrayTo_ijk(d1_gamma23_lll,ijk);
        dGamma_lll[{2,1,1}] = InterpolateArrayTo_ijk(d2_gamma11_lll,ijk);
        dGamma_lll[{2,2,2}] = InterpolateArrayTo_ijk(d2_gamma22_lll,ijk);
        dGamma_lll[{2,3,3}] = InterpolateArrayTo_ijk(d2_gamma33_lll,ijk);
        dGamma_lll[{2,1,2}] = dGamma_lll[{2,2,1}] = InterpolateArrayTo_ijk(d2_gamma12_lll,ijk);
        dGamma_lll[{2,1,3}] = dGamma_lll[{2,3,1}] = InterpolateArrayTo_ijk(d2_gamma13_lll,ijk);
        dGamma_lll[{2,2,3}] = dGamma_lll[{2,3,2}] = InterpolateArrayTo_ijk(d2_gamma23_lll,ijk);
        dGamma_lll[{3,1,1}] = InterpolateArrayTo_ijk(d3_gamma11_lll,ijk);
        dGamma_lll[{3,2,2}] = InterpolateArrayTo_ijk(d3_gamma22_lll,ijk);
        dGamma_lll[{3,3,3}] = InterpolateArrayTo_ijk(d3_gamma33_lll,ijk);
        dGamma_lll[{3,1,2}] = dGamma_lll[{3,2,1}] = InterpolateArrayTo_ijk(d3_gamma12_lll,ijk);
        dGamma_lll[{3,1,3}] = dGamma_lll[{3,3,1}] = InterpolateArrayTo_ijk(d3_gamma13_lll,ijk);
        dGamma_lll[{3,2,3}] = dGamma_lll[{3,3,2}] = InterpolateArrayTo_ijk(d3_gamma23_lll,ijk);
        return dGamma_lll;
    }
}

Tensor3x3x3 Metric::GetDerivGamma_luu(int ijk)
{
    return Tensor3x3x3
    (d1_gamma11_luu[ijk], d1_gamma12_luu[ijk], d1_gamma13_luu[ijk],
     d1_gamma12_luu[ijk], d1_gamma22_luu[ijk], d1_gamma23_luu[ijk],
     d1_gamma13_luu[ijk], d1_gamma23_luu[ijk], d1_gamma33_luu[ijk],
     d2_gamma11_luu[ijk], d2_gamma12_luu[ijk], d2_gamma13_luu[ijk],
     d2_gamma12_luu[ijk], d2_gamma22_luu[ijk], d2_gamma23_luu[ijk],
     d2_gamma13_luu[ijk], d2_gamma23_luu[ijk], d2_gamma33_luu[ijk],
     d3_gamma11_luu[ijk], d3_gamma12_luu[ijk], d3_gamma13_luu[ijk],
     d3_gamma12_luu[ijk], d3_gamma22_luu[ijk], d3_gamma23_luu[ijk],
     d3_gamma13_luu[ijk], d3_gamma23_luu[ijk], d3_gamma33_luu[ijk]);
}
Tensor3x3x3 Metric::GetDerivGamma_luu(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3x3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor3x3x3 dGamma_luu;
        dGamma_luu[{1,1,1}] = InterpolateArrayTo_ijk(d1_gamma11_luu,ijk);
        dGamma_luu[{1,2,2}] = InterpolateArrayTo_ijk(d1_gamma22_luu,ijk);
        dGamma_luu[{1,3,3}] = InterpolateArrayTo_ijk(d1_gamma33_luu,ijk);
        dGamma_luu[{1,1,2}] = dGamma_luu[{1,2,1}] = InterpolateArrayTo_ijk(d1_gamma12_luu,ijk);
        dGamma_luu[{1,1,3}] = dGamma_luu[{1,3,1}] = InterpolateArrayTo_ijk(d1_gamma13_luu,ijk);
        dGamma_luu[{1,2,3}] = dGamma_luu[{1,3,2}] = InterpolateArrayTo_ijk(d1_gamma23_luu,ijk);
        dGamma_luu[{2,1,1}] = InterpolateArrayTo_ijk(d2_gamma11_luu,ijk);
        dGamma_luu[{2,2,2}] = InterpolateArrayTo_ijk(d2_gamma22_luu,ijk);
        dGamma_luu[{2,3,3}] = InterpolateArrayTo_ijk(d2_gamma33_luu,ijk);
        dGamma_luu[{2,1,2}] = dGamma_luu[{2,2,1}] = InterpolateArrayTo_ijk(d2_gamma12_luu,ijk);
        dGamma_luu[{2,1,3}] = dGamma_luu[{2,3,1}] = InterpolateArrayTo_ijk(d2_gamma13_luu,ijk);
        dGamma_luu[{2,2,3}] = dGamma_luu[{2,3,2}] = InterpolateArrayTo_ijk(d2_gamma23_luu,ijk);
        dGamma_luu[{3,1,1}] = InterpolateArrayTo_ijk(d3_gamma11_luu,ijk);
        dGamma_luu[{3,2,2}] = InterpolateArrayTo_ijk(d3_gamma22_luu,ijk);
        dGamma_luu[{3,3,3}] = InterpolateArrayTo_ijk(d3_gamma33_luu,ijk);
        dGamma_luu[{3,1,2}] = dGamma_luu[{3,2,1}] = InterpolateArrayTo_ijk(d3_gamma12_luu,ijk);
        dGamma_luu[{3,1,3}] = dGamma_luu[{3,3,1}] = InterpolateArrayTo_ijk(d3_gamma13_luu,ijk);
        dGamma_luu[{3,2,3}] = dGamma_luu[{3,3,2}] = InterpolateArrayTo_ijk(d3_gamma23_luu,ijk);
        return dGamma_luu;
    }
}

Tensor3x3 Metric::GetExtrCurv_ll(int ijk)
{
    return Tensor3x3
    (K11_ll[ijk],K12_ll[ijk],K13_ll[ijk],
     K12_ll[ijk],K22_ll[ijk],K23_ll[ijk],
     K13_ll[ijk],K23_ll[ijk],K33_ll[ijk]);
}
Tensor3x3 Metric::GetExtrCurv_ll(const Coord& xyz)
{
    if(grid.OutsideDomain(xyz))
        return Tensor3x3(0.0);
    else
    {
        Coord ijk = grid.ijk(xyz);
        Tensor3x3 K_ll;
        K_ll[{1,1}] = InterpolateArrayTo_ijk(K11_ll,ijk);
        K_ll[{2,2}] = InterpolateArrayTo_ijk(K22_ll,ijk);
        K_ll[{3,3}] = InterpolateArrayTo_ijk(K33_ll,ijk);
        K_ll[{1,2}] = K_ll[{2,1}] = InterpolateArrayTo_ijk(K12_ll,ijk);
        K_ll[{1,3}] = K_ll[{3,1}] = InterpolateArrayTo_ijk(K13_ll,ijk);
        K_ll[{2,3}] = K_ll[{3,2}] = InterpolateArrayTo_ijk(K23_ll,ijk);
        return K_ll;
    }
}