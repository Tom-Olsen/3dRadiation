#include <iostream>
#include <span>
#include <vector>
#include "../src/Radiation.h"

using namespace std;

Tensor4x4 GrammSchmidt(const Tensor4x4 &g_ll, const Tensor4 &v0, const Tensor4 &v1, const Tensor4 &v2, const Tensor4 &v3)
{
    Tensor4 u0 = v0;
    double u0u0 = Dot(u0, u0, g_ll);
    Tensor4 e0 = u0 / sqrt(abs(u0u0));

    Tensor4 u1 = v1 - u0 * Dot(v1, u0, g_ll) / u0u0;
    double u1u1 = Dot(u1, u1, g_ll);
    Tensor4 e1 = u1 / sqrt(abs(u1u1));

    Tensor4 u2 = v2 - u0 * Dot(v2, u0, g_ll) / u0u0 - u1 * Dot(v2, u1, g_ll) / u1u1;
    double u2u2 = Dot(u2, u2, g_ll);
    Tensor4 e2 = u2 / sqrt(abs(u2u2));

    Tensor4 u3 = v3 - u0 * Dot(v3, u0, g_ll) / u0u0 - u1 * Dot(v3, u1, g_ll) / u1u1 - u2 * Dot(v3, u2, g_ll) / u2u2;
    double u3u3 = Dot(u3, u3, g_ll);
    Tensor4 e3 = u3 / sqrt(abs(u3u3));

    return Tensor4x4(e0[0], e1[0], e2[0], e3[0],
                     e0[1], e1[1], e2[1], e3[1],
                     e0[2], e1[2], e2[2], e3[2],
                     e0[3], e1[3], e2[3], e3[3]);
}

Tensor4x4 MyMethod(const Tensor4x4 &g_ll, const Tensor4 &n0)
{
    // Get 2+1 components:
    double det = (g_ll[{2, 3}] * g_ll[{2, 3}] - g_ll[{2, 2}] * g_ll[{3, 3}]);
    double a = sqrt(g_ll[{1, 1}] + (g_ll[{1, 2}] * g_ll[{1, 2}] * g_ll[{3, 3}] + g_ll[{1, 3}] * g_ll[{1, 3}] * g_ll[{2, 2}] - 2 * g_ll[{1, 2}] * g_ll[{1, 3}] * g_ll[{2, 3}]) / det);
    double bx = (g_ll[{1, 2}] * g_ll[{3, 3}] - g_ll[{1, 3}] * g_ll[{2, 3}]) / det;
    double by = (g_ll[{1, 3}] * g_ll[{2, 2}] - g_ll[{1, 2}] * g_ll[{2, 3}]) / det;

    // Get 1+1 components:
    double A = sqrt(g_ll[{2, 2}] - g_ll[{2, 3}] * g_ll[{2, 3}] / g_ll[{3, 3}]);
    double B = -g_ll[{2, 3}] / g_ll[{3, 3}];

    // Build eulerian observer like velocities:
    Tensor3 n1(1.0 / a, bx / a, by / a);
    Tensor2 n2(1.0 / A, B / A);
    double n3 = 1.0 / sqrt(g_ll[{3, 3}]);

    return Tensor4x4(n0[0], 0, 0, 0,
                     n0[1], n1[1], 0, 0,
                     n0[2], n1[2], n2[2], 0,
                     n0[3], n1[3], n2[3], n3);
}
Tensor4x4 MyMethodBoosted(Tensor4x4 &g_ll, const Tensor4 &n0)
{
    // Get 2+1 components:
    double det = (g_ll[{2, 3}] * g_ll[{2, 3}] - g_ll[{2, 2}] * g_ll[{3, 3}]);
    double a = sqrt(g_ll[{1, 1}] + (g_ll[{1, 2}] * g_ll[{1, 2}] * g_ll[{3, 3}] + g_ll[{1, 3}] * g_ll[{1, 3}] * g_ll[{2, 2}] - 2 * g_ll[{1, 2}] * g_ll[{1, 3}] * g_ll[{2, 3}]) / det);
    double bx = (g_ll[{1, 2}] * g_ll[{3, 3}] - g_ll[{1, 3}] * g_ll[{2, 3}]) / det;
    double by = (g_ll[{1, 3}] * g_ll[{2, 2}] - g_ll[{1, 2}] * g_ll[{2, 3}]) / det;

    // Get 1+1 components:
    double A = sqrt(g_ll[{2, 2}] - g_ll[{2, 3}] * g_ll[{2, 3}] / g_ll[{3, 3}]);
    double B = -g_ll[{2, 3}] / g_ll[{3, 3}];

    // Build eulerian observer like velocities:
    Tensor3 n1(1.0 / a, bx / a, by / a);
    Tensor2 n2(1.0 / A, B / A);
    double n3 = 1.0 / sqrt(g_ll[{3, 3}]);

    Tensor4x4 tetrad(n0[0], 0, 0, 0,
                     n0[1], n1[1], 0, 0,
                     n0[2], n1[2], n2[2], 0,
                     n0[3], n1[3], n2[3], n3);

    // Lorentz-Boost:
    Tensor4x4 boost(0.0);
    double vx = n0[1] / tetrad[{1, 1}];
    double vy = -n0[1] * tetrad[{2, 1}] / (tetrad[{1, 1}] * tetrad[{2, 2}]) + n0[2] / tetrad[{2, 2}];
    double vz = n0[1] * (tetrad[{2, 1}] * tetrad[{3, 2}] - tetrad[{3, 1}] * tetrad[{2, 2}]) / (tetrad[{1, 1}] * tetrad[{2, 2}] * tetrad[{3, 3}]) - n0[2] * tetrad[{3, 2}] / (tetrad[{2, 2}] * tetrad[{3, 3}]) + n0[3] / tetrad[{3, 3}];

    Tensor4 v(0, vx, vy, vz);
    double vv = v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
    double gamma = pow(1.0 - vv, -0.5); // = 1.0 / (tetrad[{0,0}] * sqrt(-g_ll[{0,0}]));
    boost[{0, 0}] = gamma;
    if (vv == 0)
        boost = Tensor4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    else
        for (int i = 1; i < 4; i++)
        {
            boost[{0, i}] = boost[{i, 0}] = -gamma * v[i];
            for (int j = 1; j < 4; j++)
                boost[{i, j}] = (gamma - 1) * v[i] * v[j] / vv + (i == j);
        }

    Tensor4x4 boostedTetrad(0);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int J = 0; J < 4; J++)
                boostedTetrad[{i, J}] += tetrad[{i, j}] * boost[{j, J}];
    return boostedTetrad;
}

void GrammSchmidtTest()
{
    size_t nx = 20;
    size_t ny = 20;
    size_t nz = 20;
    Coord start(-2, -2, -2);
    Coord end(2, 2, 2);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Tensor4 v0(1, 0, 0, 0);
    Tensor4 v1(0, 1, 0, 0);
    Tensor4 v2(0, 0, 1, 0);
    Tensor4 v3(0, 0, 0, 1);
    size_t j = 18;
    size_t k = 18;
    for (size_t i = 0; i < nx; i++)
    {
        size_t ijk = grid.Index(i, j, k);
        Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
        Tensor4x4 tetrad = GrammSchmidt(g_ll, v0, v1, v2, v3);
        cout << "i: " << i << endl;
        tetrad.Print("tetrad");
        // Tensor4x4 eta(0);
        // for (int a = 0; a < 4; a++)
        // for (int b = 0; b < 4; b++)
        // for (int A = 0; A < 4; A++)
        // for (int B = 0; B < 4; B++)
        // eta[{A, B}] += g_ll[{a, b}] * tetrad[{a, A}] * tetrad[{b, B}];
        // eta.Print("eta");
    }
}

void MyMethodTest()
{
    size_t nx = 20;
    size_t ny = 20;
    size_t nz = 20;
    Coord start(-2, -2, -2);
    Coord end(2, 2, 2);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    size_t j = 18;
    size_t k = 18;
    for (size_t i = 0; i < nx; i++)
    {
        size_t ijk = grid.Index(i, j, k);
        Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
        Tensor4 n0 = metric.uEulObs(ijk);
        // Tensor4x4 tetrad = MyMethod(g_ll, n0);
        Tensor4x4 tetrad = MyMethodBoosted(g_ll, n0);
        cout << "i: " << i << endl;
        tetrad.Print("tetrad");
        // Tensor4x4 eta(0);
        // for (int a = 0; a < 4; a++)
        // for (int b = 0; b < 4; b++)
        // for (int A = 0; A < 4; A++)
        // for (int B = 0; B < 4; B++)
        // eta[{A, B}] += g_ll[{a, b}] * tetrad[{a, A}] * tetrad[{b, B}];
        // eta.Print("eta");
    }
}

void TetradBenchmark()
{
    size_t nx = 200;
    size_t ny = 200;
    size_t nz = 200;
    Coord start(2, 2, 2);
    Coord end(3, 3, 3);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.5);
    int iterations = 1000;

    Profiler::Session &session = Profiler::Session::Get();
    session.Start("config.name", OUTPUTDIR + "testProfileResults.json");
    // Gramm Schmidt Method:
    Tensor4x4 tetrad;
    Tensor4 v0(1, 0, 0, 0);
    Tensor4 v1(0, 1, 0, 0);
    Tensor4 v2(0, 0, 1, 0);
    Tensor4 v3(0, 0, 0, 1);
    for (int i = 0; i < iterations; i++)
    {
        PROFILE_SCOPE("GS stationary");
        for (int ijk = 0; ijk < grid.nxyz; ijk++)
        {
            Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
            tetrad = GrammSchmidt(g_ll, v0, v1, v2, v3);
        }
    }
    for (int i = 0; i < iterations; i++)
    {
        PROFILE_SCOPE("GS comoving");
        for (int ijk = 0; ijk < grid.nxyz; ijk++)
        {
            Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
            Tensor4 n = metric.uEulObs(ijk);
            tetrad = GrammSchmidt(g_ll, n, v1, v2, v3);
        }
    }

    // My Methdod:
    for (int i = 0; i < iterations; i++)
    {
        PROFILE_SCOPE("Sep comoving");
        for (int ijk = 0; ijk < grid.nxyz; ijk++)
        {
            Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
            Tensor4 n0 = metric.uEulObs(ijk);
            tetrad = MyMethod(g_ll, n0);
        }
    }
    for (int i = 0; i < iterations; i++)
    {
        PROFILE_SCOPE("Sep stationary");
        for (int ijk = 0; ijk < grid.nxyz; ijk++)
        {
            Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
            Tensor4 n0 = metric.uEulObs(ijk);
            tetrad = MyMethodBoosted(g_ll, n0);
        }
    }
    session.End();

    std::vector<std::string> names = session.GetAllFunctionNames();
    for (int i = 0; i < names.size(); i++)
        session.PrintFunctionDuration(names[i]);
}

void InterpolationBenchmark()
{
    size_t nx = 200;
    size_t ny = 200;
    size_t nz = 200;
    Coord start(-2, 2.1, -2);
    Coord end(2, 3, 2);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Tensor4x4 *tetrads0 = new Tensor4x4[nx * ny * nz]();
    Tensor4x4 *tetrads1 = new Tensor4x4[nx * ny * nz]();
    Tensor4x4 *tetrads2 = new Tensor4x4[nx * ny * nz]();
    Tensor4x4 *tetrads3 = new Tensor4x4[nx * ny * nz]();

    Tensor4 v0(1, 0, 0, 0);
    Tensor4 v1(0, 1, 0, 0);
    Tensor4 v2(0, 0, 1, 0);
    Tensor4 v3(0, 0, 0, 1);
    for (size_t k = 0; k < nz; k++)
        for (size_t j = 0; j < ny; j++)
            for (size_t i = 0; i < nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Tensor4x4 g_ll = metric.GetMetric_ll(ijk);
                Tensor4 n0 = metric.uEulObs(ijk);
                tetrads0[ijk] = MyMethod(g_ll, n0);
                tetrads1[ijk] = MyMethodBoosted(g_ll, n0);
                tetrads2[ijk] = GrammSchmidt(g_ll, v0, v1, v2, v3);
                tetrads3[ijk] = GrammSchmidt(g_ll, metric.uEulObs(ijk), v1, v2, v3);
            }

    double maxNorm0 = 0;
    double maxNorm1 = 0;
    double maxNorm2 = 0;
    double maxNorm3 = 0;
    for (double k = 0.5; k < nz - 1; k++)
        for (double j = 0.5; j < ny - 1; j++)
            for (double i = 0.5; i < nx - 1; i++)
            {
                Coord xyz = grid.xyz(i, j, k);
                int i0 = floor(i);
                int i1 = ceil(i);
                int j0 = floor(j);
                int j1 = ceil(j);
                int k0 = floor(k);
                int k1 = ceil(k);
                Tensor4x4 tetrad0_000 = tetrads0[grid.Index(i0, j0, k0)];
                Tensor4x4 tetrad1_000 = tetrads1[grid.Index(i0, j0, k0)];
                Tensor4x4 tetrad0_001 = tetrads0[grid.Index(i0, j0, k1)];
                Tensor4x4 tetrad1_001 = tetrads1[grid.Index(i0, j0, k1)];
                Tensor4x4 tetrad0_010 = tetrads0[grid.Index(i0, j1, k0)];
                Tensor4x4 tetrad1_010 = tetrads1[grid.Index(i0, j1, k0)];
                Tensor4x4 tetrad0_011 = tetrads0[grid.Index(i0, j1, k1)];
                Tensor4x4 tetrad1_011 = tetrads1[grid.Index(i0, j1, k1)];
                Tensor4x4 tetrad0_100 = tetrads0[grid.Index(i1, j0, k0)];
                Tensor4x4 tetrad1_100 = tetrads1[grid.Index(i1, j0, k0)];
                Tensor4x4 tetrad0_101 = tetrads0[grid.Index(i1, j0, k1)];
                Tensor4x4 tetrad1_101 = tetrads1[grid.Index(i1, j0, k1)];
                Tensor4x4 tetrad0_110 = tetrads0[grid.Index(i1, j1, k0)];
                Tensor4x4 tetrad1_110 = tetrads1[grid.Index(i1, j1, k0)];
                Tensor4x4 tetrad0_111 = tetrads0[grid.Index(i1, j1, k1)];
                Tensor4x4 tetrad1_111 = tetrads1[grid.Index(i1, j1, k1)];

                Tensor4x4 tetrad2_000 = tetrads2[grid.Index(i0, j0, k0)];
                Tensor4x4 tetrad3_000 = tetrads3[grid.Index(i0, j0, k0)];
                Tensor4x4 tetrad2_001 = tetrads2[grid.Index(i0, j0, k1)];
                Tensor4x4 tetrad3_001 = tetrads3[grid.Index(i0, j0, k1)];
                Tensor4x4 tetrad2_010 = tetrads2[grid.Index(i0, j1, k0)];
                Tensor4x4 tetrad3_010 = tetrads3[grid.Index(i0, j1, k0)];
                Tensor4x4 tetrad2_011 = tetrads2[grid.Index(i0, j1, k1)];
                Tensor4x4 tetrad3_011 = tetrads3[grid.Index(i0, j1, k1)];
                Tensor4x4 tetrad2_100 = tetrads2[grid.Index(i1, j0, k0)];
                Tensor4x4 tetrad3_100 = tetrads3[grid.Index(i1, j0, k0)];
                Tensor4x4 tetrad2_101 = tetrads2[grid.Index(i1, j0, k1)];
                Tensor4x4 tetrad3_101 = tetrads3[grid.Index(i1, j0, k1)];
                Tensor4x4 tetrad2_110 = tetrads2[grid.Index(i1, j1, k0)];
                Tensor4x4 tetrad3_110 = tetrads3[grid.Index(i1, j1, k0)];
                Tensor4x4 tetrad2_111 = tetrads2[grid.Index(i1, j1, k1)];
                Tensor4x4 tetrad3_111 = tetrads3[grid.Index(i1, j1, k1)];

                Tensor4x4 tetrad0;
                Tensor4x4 tetrad1;
                Tensor4x4 tetrad2;
                Tensor4x4 tetrad3;
                for (int a = 0; a < 4; a++)
                    for (int b = 0; b < 4; b++)
                    {
                        tetrad0[{a, b}] = TrilinearInterpolation(i - i0, j - j0, k - k0,
                                                                 tetrad0_000[{a, b}], tetrad0_001[{a, b}], tetrad0_010[{a, b}], tetrad0_011[{a, b}],
                                                                 tetrad0_100[{a, b}], tetrad0_101[{a, b}], tetrad0_110[{a, b}], tetrad0_111[{a, b}]);
                        tetrad1[{a, b}] = TrilinearInterpolation(i - i0, j - j0, k - k0,
                                                                 tetrad1_000[{a, b}], tetrad1_001[{a, b}], tetrad1_010[{a, b}], tetrad1_011[{a, b}],
                                                                 tetrad1_100[{a, b}], tetrad1_101[{a, b}], tetrad1_110[{a, b}], tetrad1_111[{a, b}]);
                        tetrad2[{a, b}] = TrilinearInterpolation(i - i0, j - j0, k - k0,
                                                                 tetrad2_000[{a, b}], tetrad2_001[{a, b}], tetrad2_010[{a, b}], tetrad2_011[{a, b}],
                                                                 tetrad2_100[{a, b}], tetrad2_101[{a, b}], tetrad2_110[{a, b}], tetrad2_111[{a, b}]);
                        tetrad3[{a, b}] = TrilinearInterpolation(i - i0, j - j0, k - k0,
                                                                 tetrad3_000[{a, b}], tetrad3_001[{a, b}], tetrad3_010[{a, b}], tetrad3_011[{a, b}],
                                                                 tetrad3_100[{a, b}], tetrad3_101[{a, b}], tetrad3_110[{a, b}], tetrad3_111[{a, b}]);
                    }

                Tensor4x4 g_ll = metric.GetMetric_ll(xyz);
                Tensor4x4 etaExact = metric.GetMinkowskiMetric_ll(xyz);
                Tensor4x4 etaApprox0(0);
                Tensor4x4 etaApprox1(0);
                Tensor4x4 etaApprox2(0);
                Tensor4x4 etaApprox3(0);
                for (int a = 0; a < 4; a++)
                    for (int b = 0; b < 4; b++)
                        for (int A = 0; A < 4; A++)
                            for (int B = 0; B < 4; B++)
                            {
                                etaApprox0[{A, B}] += g_ll[{a, b}] * tetrad0[{a, A}] * tetrad0[{b, B}];
                                etaApprox1[{A, B}] += g_ll[{a, b}] * tetrad1[{a, A}] * tetrad1[{b, B}];
                                etaApprox2[{A, B}] += g_ll[{a, b}] * tetrad2[{a, A}] * tetrad2[{b, B}];
                                etaApprox3[{A, B}] += g_ll[{a, b}] * tetrad3[{a, A}] * tetrad3[{b, B}];
                            }

                double norm0 = 0;
                double norm1 = 0;
                double norm2 = 0;
                double norm3 = 0;
                for (int a = 0; a < 4; a++)
                    for (int b = 0; b < 4; b++)
                    {
                        norm0 += IntegerPow<2>(etaApprox0[{a, b}] - etaExact[{a, b}]);
                        norm1 += IntegerPow<2>(etaApprox1[{a, b}] - etaExact[{a, b}]);
                        norm2 += IntegerPow<2>(etaApprox2[{a, b}] - etaExact[{a, b}]);
                        norm3 += IntegerPow<2>(etaApprox3[{a, b}] - etaExact[{a, b}]);
                    }
                norm0 = sqrt(norm0);
                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                norm3 = sqrt(norm3);
                maxNorm0 = max(maxNorm0, norm0);
                maxNorm1 = max(maxNorm1, norm1);
                maxNorm2 = max(maxNorm2, norm2);
                maxNorm3 = max(maxNorm3, norm3);
            }
    cout << "Separation Comoving:      " << Format(maxNorm0, 12) << endl;
    cout << "Separation Stationary:    " << Format(maxNorm1, 12) << endl;
    cout << "Gramm Schmidt Stationary: " << Format(maxNorm2, 12) << endl;
    cout << "Gramm Schmidt Comoving:   " << Format(maxNorm3, 12) << endl;
}

void RealBufferTest()
{
    RealBuffer a;
    a.resize(3);
    a[0] = 0;
    a[1] = 1;
    a[2] = 2;
    a.resize(5);
    a.push_back(3);
    for (int i = 0; i < a.size(); i++)
        cout << a[i] << endl;
}

void TestInitialDataFunction()
{
    LebedevStencil stencil(21, 0.25, 0.00, 0.00);
    LebedevStencil stencilBase(21);

    double E = 1;
    double F = 0.0;
    double sigma = stencilBase.sigmaOfRelativeFlux.Evaluate(F / E);

    RealBuffer IBase;
    IBase.resize(stencilBase.nDir);
    for (int d = 0; d < stencilBase.nDir; d++)
    {
        IBase[d] = Intensity(sigma, E, stencilBase.Theta(d)) / (stencilBase.W(d) * stencilBase.nDir);
        // std::cout << IBase[d] * stencilBase.Cx(d) << "," << IBase[d] * stencilBase.Cy(d) << "," << IBase[d] * stencilBase.Cz(d) << std::endl;
    }

    RealBuffer I;
    I.resize(stencil.nDir);
    for( int d = 0; d < stencil.nDir; d++)
    {
        Vector3 dir = stencil.Cv3(d);
        std::tuple<std::vector<size_t>, std::vector<double>> neighboursAndWeights = stencilBase.VoronoiNeighboursAndWeights(dir);
        std::span<const size_t> neighbours = std::get<0>(neighboursAndWeights);
        std::span<const double> weights = std::get<1>(neighboursAndWeights);
        double interpolatetValue = 0;
        for (size_t p = 0; p < weights.size(); p++)
            interpolatetValue += weights[p] * IBase[neighbours[p]];
        I[d] = interpolatetValue;
        std::cout << I[d] * stencil.Cx(d) << "," << I[d] * stencil.Cy(d) << "," << I[d] * stencil.Cz(d) << std::endl;
    }
}



int main()
{
    // GrammSchmidtTest();
    // MyMethodTest();
    // MyMethodBoosted();
    // TetradBenchmark();
    // InterpolationBenchmark();
    TestInitialDataFunction();

    // LebedevStencil stencilA(21);
    // cout << "Lebedev21" << endl;
    // cout << "sigmaMax: " << stencilA.sigmaMax << std::endl;
    // cout << "relativeFluxMax: " << stencilA.relativeFluxMax << std::endl;

    // LebedevStencil stencilB(21,0.13, 0.0, 0.0);
    // cout << "Lebedev21 0.13 0.00 0.00" << endl;
    // cout << "sigmaMax: " << stencilB.sigmaMax << std::endl;
    // cout << "relativeFluxMax: " << stencilB.relativeFluxMax << std::endl;
}