#include "Radiation.h"

Radiation::Radiation(Metric &metric, Stencil &stencil, LebedevStencil &streamingStencil, InterpolationGrid &interpGrid, Camera &camera, Config config)
    : grid(metric.grid), metric(metric), stencil(stencil), streamingStencil(streamingStencil), interpGrid(interpGrid), camera(camera), config(config), logger(stencil, streamingStencil, metric)
{
    isInitialGridPoint = new bool[grid.nxyz]();
    initialE_LF.resize(grid.nxyz);
    initialFx_LF.resize(grid.nxyz);
    initialFy_LF.resize(grid.nxyz);
    initialFz_LF.resize(grid.nxyz);
    initialPxx_LF.resize(grid.nxyz);
    initialPxy_LF.resize(grid.nxyz);
    initialPxz_LF.resize(grid.nxyz);
    initialPyy_LF.resize(grid.nxyz);
    initialPyz_LF.resize(grid.nxyz);
    initialPzz_LF.resize(grid.nxyz);
    initialKappa0.resize(grid.nxyz);
    initialKappa1.resize(grid.nxyz);
    initialKappaA.resize(grid.nxyz);
    initialEta.resize(grid.nxyz);
    initialI.resize(grid.nxyz * stencil.nDir);
    initialQ.resize(grid.nxyz);

    q.resize(grid.nxyz);
    qNew.resize(grid.nxyz);
    E.resize(grid.nxyz);
    Fx.resize(grid.nxyz);
    Fy.resize(grid.nxyz);
    Fz.resize(grid.nxyz);
    Pxx.resize(grid.nxyz);
    Pxy.resize(grid.nxyz);
    Pxz.resize(grid.nxyz);
    Pyy.resize(grid.nxyz);
    Pyz.resize(grid.nxyz);
    Pzz.resize(grid.nxyz);
    E_LF.resize(grid.nxyz);
    Fx_LF.resize(grid.nxyz);
    Fy_LF.resize(grid.nxyz);
    Fz_LF.resize(grid.nxyz);
    Pxx_LF.resize(grid.nxyz);
    Pxy_LF.resize(grid.nxyz);
    Pxz_LF.resize(grid.nxyz);
    Pyy_LF.resize(grid.nxyz);
    Pyz_LF.resize(grid.nxyz);
    Pzz_LF.resize(grid.nxyz);
    kappa0.resize(grid.nxyz);
    kappa1.resize(grid.nxyz);
    kappaA.resize(grid.nxyz);
    eta.resize(grid.nxyz);
    I.resize(grid.nxyz * stencil.nDir);
    Inew.resize(grid.nxyz * stencil.nDir);
    coefficientsS.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsX.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsY.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsZ.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsCx.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsCy.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsCz.resize(grid.nxyz * streamingStencil.nCoefficients);

    // Initialize all Quaternions to identity:
    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
        q[ijk] = qNew[ijk] = initialQ[ijk] = glm::quat(1, 0, 0, 0);
}
Radiation::~Radiation()
{
    delete[] isInitialGridPoint;
}

size_t Radiation::Index(size_t ijk, size_t d)
{
    return d + ijk * stencil.nDir;
}
size_t Radiation::Index(size_t i, size_t j, size_t k, size_t d)
{
    size_t ijk = grid.Index(i, j, k);
    return d + ijk * stencil.nDir;
}
size_t Radiation::HarmonicIndex(size_t f, size_t ijk)
{
    return f + ijk * streamingStencil.nCoefficients;
}

Tensor4 Radiation::InitialDataLFtoIF(size_t ijk)
{
    // Transform given initial E_LF, Fx_LF, Fy_LF into initial E_IF, Fx_IF, Fy_IF:
    Tensor4x4 EnergyMomentumTensorLF(initialE_LF[ijk], initialFx_LF[ijk], initialFy_LF[ijk], initialFz_LF[ijk],
                                     initialFx_LF[ijk], initialPxx_LF[ijk], initialPxy_LF[ijk], initialPxz_LF[ijk],
                                     initialFy_LF[ijk], initialPxy_LF[ijk], initialPyy_LF[ijk], initialPyz_LF[ijk],
                                     initialFz_LF[ijk], initialPxz_LF[ijk], initialPyz_LF[ijk], initialPzz_LF[ijk]);
    Tensor4x4 EnergyMomentumTensorIF = TransformLFtoIF(EnergyMomentumTensorLF, metric.GetTetrad(ijk).Invert());

    double initialE_IF = EnergyMomentumTensorIF[{0, 0}];
    double initialFx_IF = EnergyMomentumTensorLF[{0, 1}];
    double initialFy_IF = EnergyMomentumTensorLF[{0, 2}];
    double initialFz_IF = EnergyMomentumTensorLF[{0, 3}];
    double initialPxx_IF = EnergyMomentumTensorLF[{1, 1}];
    double initialPxy_IF = EnergyMomentumTensorLF[{1, 2}];
    double initialPxz_IF = EnergyMomentumTensorLF[{1, 3}];
    double initialPyy_IF = EnergyMomentumTensorLF[{2, 2}];
    double initialPyz_IF = EnergyMomentumTensorLF[{2, 3}];
    double initialPzz_IF = EnergyMomentumTensorLF[{3, 3}];

    // |F| > E is unphysical:
    double initialF_IF = Tensor3(initialFx_IF, initialFy_IF, initialFz_IF).EuklNorm();
    if (initialF_IF > initialE_IF)
    {
        initialFx_IF *= initialE_IF / initialF_IF;
        initialFy_IF *= initialE_IF / initialF_IF;
        initialFz_IF *= initialE_IF / initialF_IF;
    }

    return Tensor4(initialE_IF, initialFx_IF, initialFy_IF, initialFz_IF);
}

void Radiation::LoadInitialData()
{
    PROFILE_FUNCTION();
    bool isAdaptiveStreaming = (config.streamingType == StreamingType::FlatAdaptive || config.streamingType == StreamingType::CurvedAdaptive);

    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
    {
        if (!isInitialGridPoint[ijk])
            continue;

        // Convert given LF initial data to IF:
        Tensor4 initialDataIF = InitialDataLFtoIF(ijk);
        double initialE_IF = initialDataIF[0];
        double initialFx_IF = initialDataIF[1];
        double initialFy_IF = initialDataIF[2];
        double initialFz_IF = initialDataIF[3];

        // Flux direction and magnitude in IF:
        Tensor3 initialFxyz_IF(initialFx_IF, initialFy_IF, initialFz_IF);
        double initialF_IF = initialFxyz_IF.EuklNorm();
        double relativeF_IF = initialF_IF / initialE_IF;
        Tensor3 dirInitialF = (isAdaptiveStreaming) ? Tensor3(0, 0, 1) : Tensor3(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF, initialFz_IF / initialF_IF);
        if (initialF_IF < MIN_FLUX_NORM)
            dirInitialF = Tensor3(0, 0, 1);

        kappa0[ijk] = initialKappa0[ijk];
        kappa1[ijk] = initialKappa1[ijk];
        kappaA[ijk] = initialKappaA[ijk];
        eta[ijk] = initialEta[ijk];
        double sigma = stencil.fluxToSigmaTable.Evaluate(relativeF_IF);
        double normalization = stencil.fluxToNormalizationTable.Evaluate(relativeF_IF);

        if (config.initialDataType == InitialDataType::Moments)
        {
            // Von Mises distribution:
            // https://en.wikipedia.org/wiki/Von_Mises_distribution
            if (isInitialGridPoint[ijk])
            {
                q[ijk] = (isAdaptiveStreaming) ? initialQ[ijk] : glm::quat(1, 0, 0, 0);
                for (size_t d = 0; d < stencil.nDir; d++)
                    I[Index(ijk, d)] = initialE_IF * exp(sigma * Tensor3::Dot(dirInitialF, stencil.Ct3(d)) - normalization);
            }
        }
        else if (config.initialDataType == InitialDataType::Intensities)
        {
            if (isInitialGridPoint[ijk])
            {
                q[ijk] = (isAdaptiveStreaming) ? initialQ[ijk] : glm::quat(1, 0, 0, 0);
                for (size_t d = 0; d < stencil.nDir; d++)
                    I[Index(ijk, d)] = initialE_LF[ijk] * initialI[Index(ijk, d)];
            }
        }
    }
}

void Radiation::UpdateSphericalHarmonicsCoefficients()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                double dataS[streamingStencil.nDir];
                double dataX[streamingStencil.nDir];
                double dataY[streamingStencil.nDir];
                double dataZ[streamingStencil.nDir];
                double dataCx[streamingStencil.nDir];
                double dataCy[streamingStencil.nDir];
                double dataCz[streamingStencil.nDir];
                Coord xyz0 = grid.xyz(i, j, k);
                double alpha = metric.GetAlpha(ijk);

                for (size_t d = 0; d < streamingStencil.nDir; d++)
                {
                    // Initial data for geodesic equation:
                    double s = 1;
                    Coord xyz = xyz0;
                    Tensor3 c = streamingStencil.Ct3(d);
                    Tensor4 uIF(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
                    Tensor3 vLF = Vec3ObservedByEulObs<IF, LF>(uIF, xyz, metric);

                    // Solve geodesic equation backwards:
                    if (!metric.InsideBH(xyz))
                        s *= RK45_GeodesicEquation<-1>(grid.dt, xyz, vLF, metric);
                    else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
                        vLF = Tensor3(0.0);

                    Tensor3 vIF = TransformLFtoIF(vLF, metric.GetTetradInverse(xyz));

                    // Final data points for fourier expansion:
                    dataS[d] = 1.0 / s;
                    dataX[d] = xyz[1];
                    dataY[d] = xyz[2];
                    dataZ[d] = xyz[3];
                    dataCx[d] = vIF[1];
                    dataCy[d] = vIF[2];
                    dataCz[d] = vIF[3];
                }
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataS, &coefficientsS[HarmonicIndex(0, ijk)]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataX, &coefficientsX[HarmonicIndex(0, ijk)]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataY, &coefficientsY[HarmonicIndex(0, ijk)]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataZ, &coefficientsZ[HarmonicIndex(0, ijk)]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[HarmonicIndex(0, ijk)]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[HarmonicIndex(0, ijk)]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCz, &coefficientsCz[HarmonicIndex(0, ijk)]);
            }
}

void Radiation::ComputeMomentsIF()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
    {
        E[ijk] = 0.0;
        Fx[ijk] = 0.0;
        Fy[ijk] = 0.0;
        Fz[ijk] = 0.0;
        Pxx[ijk] = 0.0;
        Pxy[ijk] = 0.0;
        Pxz[ijk] = 0.0;
        Pyy[ijk] = 0.0;
        Pyz[ijk] = 0.0;
        Pzz[ijk] = 0.0;
        for (size_t d = 0; d < stencil.nDir; d++)
        {
            // Skip ghost directions:
            if (stencil.W(d) == 0.0)
                continue;

            Tensor3 dir = q[ijk] * stencil.Ct3(d);
            size_t index = Index(ijk, d);
            double c = stencil.W(d) * I[index];

            E[ijk] += c;
            Fx[ijk] += c * dir[1];
            Fy[ijk] += c * dir[2];
            Fz[ijk] += c * dir[3];
            Pxx[ijk] += c * dir[1] * dir[1];
            Pxy[ijk] += c * dir[1] * dir[2];
            Pxz[ijk] += c * dir[1] * dir[3];
            Pyy[ijk] += c * dir[2] * dir[2];
            Pyz[ijk] += c * dir[2] * dir[3];
            Pzz[ijk] += c * dir[3] * dir[3];
        }
    }
}
void Radiation::ComputeMomentsLF()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                int ijk = grid.Index(i, j, k);
                if (metric.InsideBH(grid.xyz(i, j, k)))
                {
                    E_LF[ijk] = 0.0;
                    Fx_LF[ijk] = 0.0;
                    Fy_LF[ijk] = 0.0;
                    Fz_LF[ijk] = 0.0;
                    Pxx_LF[ijk] = 0.0;
                    Pxy_LF[ijk] = 0.0;
                    Pxz_LF[ijk] = 0.0;
                    Pyy_LF[ijk] = 0.0;
                    Pyz_LF[ijk] = 0.0;
                    Pzz_LF[ijk] = 0.0;
                    continue;
                }
                Tensor4x4 EnergyMomentumTensorIF(E[ijk], Fx[ijk], Fy[ijk], Fz[ijk],
                                                 Fx[ijk], Pxx[ijk], Pxy[ijk], Pxz[ijk],
                                                 Fy[ijk], Pxy[ijk], Pyy[ijk], Pyz[ijk],
                                                 Fz[ijk], Pxz[ijk], Pyz[ijk], Pzz[ijk]);
                Tensor4x4 EnergyMomentumTensorLF = TransformIFtoLF(EnergyMomentumTensorIF, metric.GetTetrad(ijk));

                E_LF[ijk] = EnergyMomentumTensorLF[{0, 0}];
                Fx_LF[ijk] = EnergyMomentumTensorLF[{0, 1}];
                Fy_LF[ijk] = EnergyMomentumTensorLF[{0, 2}];
                Fz_LF[ijk] = EnergyMomentumTensorLF[{0, 3}];
                Pxx_LF[ijk] = EnergyMomentumTensorLF[{1, 1}];
                Pxy_LF[ijk] = EnergyMomentumTensorLF[{1, 2}];
                Pxz_LF[ijk] = EnergyMomentumTensorLF[{1, 3}];
                Pyy_LF[ijk] = EnergyMomentumTensorLF[{2, 2}];
                Pyz_LF[ijk] = EnergyMomentumTensorLF[{2, 3}];
                Pzz_LF[ijk] = EnergyMomentumTensorLF[{3, 3}];
            }
}

Coord Radiation::GetTempCoordinate(size_t ijk, const Tensor3 &direction)
{
    Coord xyzTemp;
    xyzTemp[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsX[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
    xyzTemp[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsY[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
    xyzTemp[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsZ[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
    return xyzTemp;
}
Tensor3 Radiation::GetTemp3VelocityIF(size_t ijk, const Tensor3 &direction)
{
    Tensor3 vTempIF;
    vTempIF[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCx[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
    vTempIF[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCy[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
    vTempIF[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCz[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
    return vTempIF;
}
double Radiation::GetFrequencyShift(size_t ijk, const Tensor3 &direction)
{
    return SphericalHarmonicsXyz::GetValue(direction, &coefficientsS[HarmonicIndex(0, ijk)], streamingStencil.nCoefficients);
}

double Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF)
{
    // Interpolation point:
    vTempIF = Invert(q[ijk]) * vTempIF;

    // Barycentric interpolation (way to slow):
    // std::tuple<Vector3Int, Vector3> neighboursAndWeights = stencil.BarycentricNeighboursAndWeights(Vector3(vTempIF[1], vTempIF[2], vTempIF[3]));
    // Vector3Int neighbours = std::get<0>(neighboursAndWeights);
    // Vector3 weights = std::get<1>(neighboursAndWeights);
    // double I0 = weights[0] * I[Index(ijk, neighbours[0])];
    // double I1 = weights[1] * I[Index(ijk, neighbours[1])];
    // double I2 = weights[2] * I[Index(ijk, neighbours[2])];
    // return std::max(0.0, I0 + I1 + I2);

    // Voronoi Interpolation (way to slow):
    // std::tuple<std::vector<size_t>, std::vector<double>> neighboursAndWeights = stencil.VoronoiNeighboursAndWeights(Vector3(vTempIF[1], vTempIF[2], vTempIF[3]));
    // std::vector<size_t> neighbours = std::get<0>(neighboursAndWeights);
    // std::vector<double> weights = std::get<1>(neighboursAndWeights);
    // double value = 0.0;
    // for (int i = 0; i < neighbours.size(); i++)
    //     value += weights[i] * I[Index(ijk, neighbours[i])];
    // return std::max(0.0, value);

    // Optimized Barycentric Interpolation (weirdly slower than optimized voronoi interpolation):
    // double i = interpGrid.i(vTempIF);
    // double j = interpGrid.j(vTempIF);
    // size_t i0 = std::floor(i);
    // size_t i1 = i0 + 1;
    // size_t j0 = std::floor(j);
    // size_t j1 = (j0 + 1) % interpGrid.nPh;
    // Vector3Int neighbours00 = interpGrid.BarycentricNeighbours(i0, j0);
    // Vector3Int neighbours01 = interpGrid.BarycentricNeighbours(i0, j1);
    // Vector3Int neighbours10 = interpGrid.BarycentricNeighbours(i1, j0);
    // Vector3Int neighbours11 = interpGrid.BarycentricNeighbours(i1, j1);
    // Vector3 weights00 = interpGrid.BarycentricWeights(i0, j0);
    // Vector3 weights01 = interpGrid.BarycentricWeights(i0, j1);
    // Vector3 weights10 = interpGrid.BarycentricWeights(i1, j0);
    // Vector3 weights11 = interpGrid.BarycentricWeights(i1, j1);
    // double value00 = 0;
    // double value01 = 0;
    // double value10 = 0;
    // double value11 = 0;
    // for (int i = 0; i < 3; i++)
    // {
    //     value00 += weights00[i] * I[Index(ijk, neighbours00[i])];
    //     value01 += weights01[i] * I[Index(ijk, neighbours01[i])];
    //     value10 += weights10[i] * I[Index(ijk, neighbours10[i])];
    //     value11 += weights11[i] * I[Index(ijk, neighbours11[i])];
    // }
    // return std::max(0.0, BilinearInterpolation(i - i0, j - j0, value00, value01, value10, value11));

    // Optimized Voronoi Interpolation:
    double i = interpGrid.i(vTempIF);
    double j = interpGrid.j(vTempIF);
    size_t i0 = std::floor(i);
    size_t i1 = i0 + 1;
    size_t j0 = std::floor(j);
    size_t j1 = (j0 + 1) % interpGrid.nPh;
    std::span<const size_t> neighbours00 = interpGrid.VoronoiNeighbours(i0, j0);
    std::span<const size_t> neighbours01 = interpGrid.VoronoiNeighbours(i0, j1);
    std::span<const size_t> neighbours10 = interpGrid.VoronoiNeighbours(i1, j0);
    std::span<const size_t> neighbours11 = interpGrid.VoronoiNeighbours(i1, j1);
    std::span<const double> weights00 = interpGrid.VoronoiWeights(i0, j0);
    std::span<const double> weights01 = interpGrid.VoronoiWeights(i0, j1);
    std::span<const double> weights10 = interpGrid.VoronoiWeights(i1, j0);
    std::span<const double> weights11 = interpGrid.VoronoiWeights(i1, j1);
    double value00 = 0;
    double value01 = 0;
    double value10 = 0;
    double value11 = 0;
    for (size_t i = 0; i < weights00.size(); i++)
        value00 += weights00[i] * I[Index(ijk, neighbours00[i])];
    for (size_t i = 0; i < weights01.size(); i++)
        value01 += weights01[i] * I[Index(ijk, neighbours01[i])];
    for (size_t i = 0; i < weights10.size(); i++)
        value10 += weights10[i] * I[Index(ijk, neighbours10[i])];
    for (size_t i = 0; i < weights11.size(); i++)
        value11 += weights11[i] * I[Index(ijk, neighbours11[i])];
    return std::max(0.0, BilinearInterpolation(i - i0, j - j0, value00, value01, value10, value11));
}

Tensor3 Radiation::AverageF(size_t i, size_t j, size_t k)
{
    // i,j,k >= 2, thus i+a etc will never be negative.
    Tensor3 averageF(0.0);
    for (int c = -1; c <= 1; c++)
        for (int b = -1; b <= 1; b++)
            for (int a = -1; a <= 1; a++)
            {
                size_t index = grid.Index(i + a, j + b, k + c);
                averageF[1] += Fx[index];
                averageF[2] += Fy[index];
                averageF[3] += Fz[index];
            }
    averageF[1] /= 27.0;
    averageF[2] /= 27.0;
    averageF[3] /= 27.0;
    return averageF;
}

void Radiation::UpdateQuaternions()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Tensor3 averageF = AverageF(i, j, k);
                double norm = averageF.EuklNorm();

                // At least 1% of the light points in the direction of the first momentum
                if (E[ijk] > 1e-16 && norm / E[ijk] > 0.01)
                {
                    glm::vec3 to(averageF[1] / norm, averageF[2] / norm, averageF[3] / norm);
                    qNew[ijk] = glm::quat(from, to);
                }
                else
                    qNew[ijk] = q[ijk];
            }
}

void Radiation::StreamFlatFixed()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                // Index of lattice point ijk:
                size_t ijk = grid.Index(i, j, k);
                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    // Index of population d at lattice point ijk:
                    size_t index = Index(ijk, d);

                    // Get temp velocity:
                    Tensor3 direction = stencil.Ct3(d);

                    // Get temp lattice point:
                    Coord xyzTemp = grid.xyz(i, j, k);
                    xyzTemp[1] -= direction[1] * grid.dt;
                    xyzTemp[2] -= direction[2] * grid.dt;
                    xyzTemp[3] -= direction[3] * grid.dt;

                    // Get 8 nearest Grid Points:
                    double iTemp = grid.i(xyzTemp[1]);
                    double jTemp = grid.j(xyzTemp[2]);
                    double kTemp = grid.k(xyzTemp[3]);
                    size_t i0 = std::floor(iTemp);
                    size_t i1 = i0 + 1;
                    size_t j0 = std::floor(jTemp);
                    size_t j1 = j0 + 1;
                    size_t k0 = std::floor(kTemp);
                    size_t k1 = k0 + 1;

                    Inew[index] = TrilinearInterpolation(iTemp - i0, jTemp - j0, kTemp - k0,
                                                         I[Index(i0, j0, k0, d)], I[Index(i0, j0, k1, d)], I[Index(i0, j1, k0, d)], I[Index(i0, j1, k1, d)],
                                                         I[Index(i1, j0, k0, d)], I[Index(i1, j0, k1, d)], I[Index(i1, j1, k0, d)], I[Index(i1, j1, k1, d)]);
                }
            }
    std::swap(I, Inew);
}
void Radiation::StreamFlatAdaptive()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                // Index of lattice point ijk:
                size_t ijk = grid.Index(i, j, k);
                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    // Index of population d at lattice point ijk:
                    size_t index = Index(ijk, d);

                    // Get temp velocity:
                    Tensor3 direction = qNew[ijk] * stencil.Ct3(d);

                    // Get temp lattice point:
                    Coord xyzTemp = grid.xyz(i, j, k);
                    xyzTemp[1] -= direction[1] * grid.dt;
                    xyzTemp[2] -= direction[2] * grid.dt;
                    xyzTemp[3] -= direction[3] * grid.dt;

                    // Get 8 nearest Grid Points:
                    double iTemp = grid.i(xyzTemp[1]);
                    double jTemp = grid.j(xyzTemp[2]);
                    double kTemp = grid.k(xyzTemp[3]);
                    size_t i0 = std::floor(iTemp);
                    size_t i1 = i0 + 1;
                    size_t j0 = std::floor(jTemp);
                    size_t j1 = j0 + 1;
                    size_t k0 = std::floor(kTemp);
                    size_t k1 = k0 + 1;

                    double intensityAt_i0j0k0 = IntensityAt(grid.Index(i0, j0, k0), direction);
                    double intensityAt_i0j0k1 = IntensityAt(grid.Index(i0, j0, k1), direction);
                    double intensityAt_i0j1k0 = IntensityAt(grid.Index(i0, j1, k0), direction);
                    double intensityAt_i0j1k1 = IntensityAt(grid.Index(i0, j1, k1), direction);
                    double intensityAt_i1j0k0 = IntensityAt(grid.Index(i1, j0, k0), direction);
                    double intensityAt_i1j0k1 = IntensityAt(grid.Index(i1, j0, k1), direction);
                    double intensityAt_i1j1k0 = IntensityAt(grid.Index(i1, j1, k0), direction);
                    double intensityAt_i1j1k1 = IntensityAt(grid.Index(i1, j1, k1), direction);
                    Inew[index] = TrilinearInterpolation(iTemp - i0, jTemp - j0, kTemp - k0,
                                                         intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
                                                         intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
                }
            }
    std::swap(I, Inew);
    std::swap(q, qNew);
}
void Radiation::StreamCurvedFixed()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                // Index of lattice point ijk:
                size_t ijk = grid.Index(i, j, k);

                // Skip LPs which are inside BH:
                if (metric.InsideBH(grid.xyz(i, j, k)))
                {
                    for (size_t d = 0; d < stencil.nDir; d++)
                        Inew[Index(ijk, d)] = 0;
                    continue;
                }

                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    // Index of population d at lattice point ijk:
                    size_t index = Index(ijk, d);

                    // Get velocity direction in IF:
                    Tensor3 direction = stencil.Ct3(d);

                    // Get quantities at emission point:
                    double s = GetFrequencyShift(ijk, direction);
                    Coord xyzTemp = GetTempCoordinate(ijk, direction);
                    Tensor3 vTempIF = GetTemp3VelocityIF(ijk, direction);

                    // Skip temporary Grid Points inside BH:
                    if (metric.InsideBH(xyzTemp))
                    {
                        Inew[index] = 0;
                        continue;
                    }

                    // Get 8 nearest Grid Points:
                    double iTemp = grid.i(xyzTemp[1]);
                    double jTemp = grid.j(xyzTemp[2]);
                    double kTemp = grid.k(xyzTemp[3]);
                    size_t i0 = std::floor(iTemp);
                    size_t i1 = i0 + 1;
                    size_t j0 = std::floor(jTemp);
                    size_t j1 = j0 + 1;
                    size_t k0 = std::floor(kTemp);
                    size_t k1 = k0 + 1;

                    // Intensity interpolation:
                    double alpha = metric.GetAlpha(ijk);
                    double intensityAt_i0j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j0, k0))) * IntensityAt(grid.Index(i0, j0, k0), vTempIF);
                    double intensityAt_i0j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j0, k1))) * IntensityAt(grid.Index(i0, j0, k1), vTempIF);
                    double intensityAt_i0j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j1, k0))) * IntensityAt(grid.Index(i0, j1, k0), vTempIF);
                    double intensityAt_i0j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j1, k1))) * IntensityAt(grid.Index(i0, j1, k1), vTempIF);
                    double intensityAt_i1j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j0, k0))) * IntensityAt(grid.Index(i1, j0, k0), vTempIF);
                    double intensityAt_i1j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j0, k1))) * IntensityAt(grid.Index(i1, j0, k1), vTempIF);
                    double intensityAt_i1j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j1, k0))) * IntensityAt(grid.Index(i1, j1, k0), vTempIF);
                    double intensityAt_i1j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j1, k1))) * IntensityAt(grid.Index(i1, j1, k1), vTempIF);

                    // Interpolate intensity from neighbouring 8 lattice points to temporary point:
                    Inew[index] = IntegerPow<4>(s) * TrilinearInterpolation(iTemp - i0, jTemp - j0, kTemp - k0,
                                                                            intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
                                                                            intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
                }
            }

    std::swap(I, Inew);
}
void Radiation::StreamCurvedAdaptive()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                // Index of lattice point ijk:
                size_t ijk = grid.Index(i, j, k);

                // Skip LPs which are inside BH:
                if (metric.InsideBH(grid.xyz(i, j, k)))
                {
                    for (size_t d = 0; d < stencil.nDir; d++)
                        Inew[Index(ijk, d)] = 0;
                    continue;
                }

                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    // Index of population d at lattice point ijk:
                    size_t index = Index(ijk, d);

                    // Skip LPs which are inside BH:
                    if (metric.InsideBH(grid.xyz(i, j, k)))
                    {
                        Inew[index] = 0;
                        continue;
                    }

                    // Get velocity direction in IF:
                    Tensor3 direction = qNew[ijk] * stencil.Ct3(d);

                    // Get quantities at emission point:
                    double s = GetFrequencyShift(ijk, direction);
                    Coord xyzTemp = GetTempCoordinate(ijk, direction);
                    Tensor3 vTempIF = GetTemp3VelocityIF(ijk, direction);

                    // Skip temporary Grid Points inside BH:
                    if (metric.InsideBH(xyzTemp))
                    {
                        Inew[index] = 0;
                        continue;
                    }

                    // Get 8 nearest Grid Points:
                    double iTemp = grid.i(xyzTemp[1]);
                    double jTemp = grid.j(xyzTemp[2]);
                    double kTemp = grid.k(xyzTemp[3]);
                    size_t i0 = std::floor(iTemp);
                    size_t i1 = i0 + 1;
                    size_t j0 = std::floor(jTemp);
                    size_t j1 = j0 + 1;
                    size_t k0 = std::floor(kTemp);
                    size_t k1 = k0 + 1;

                    // Intensity interpolation:
                    double alpha = metric.GetAlpha(ijk);
                    double intensityAt_i0j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j0, k0))) * IntensityAt(grid.Index(i0, j0, k0), vTempIF);
                    double intensityAt_i0j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j0, k1))) * IntensityAt(grid.Index(i0, j0, k1), vTempIF);
                    double intensityAt_i0j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j1, k0))) * IntensityAt(grid.Index(i0, j1, k0), vTempIF);
                    double intensityAt_i0j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j1, k1))) * IntensityAt(grid.Index(i0, j1, k1), vTempIF);
                    double intensityAt_i1j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j0, k0))) * IntensityAt(grid.Index(i1, j0, k0), vTempIF);
                    double intensityAt_i1j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j0, k1))) * IntensityAt(grid.Index(i1, j0, k1), vTempIF);
                    double intensityAt_i1j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j1, k0))) * IntensityAt(grid.Index(i1, j1, k0), vTempIF);
                    double intensityAt_i1j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j1, k1))) * IntensityAt(grid.Index(i1, j1, k1), vTempIF);

                    // Interpolate intensity from neighbouring 8 lattice points to temporary point:
                    Inew[index] = IntegerPow<4>(s) * TrilinearInterpolation(iTemp - i0, jTemp - j0, kTemp - k0,
                                                                            intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
                                                                            intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
                }
            }
    std::swap(I, Inew);
    std::swap(q, qNew);
}

void Radiation::CollideStaticFluidForwardEuler()
{
    PROFILE_FUNCTION();

    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                if (metric.InsideBH(grid.xyz(i, j, k)))
                    continue;
                size_t ijk = grid.Index(i, j, k);

                double alpha = metric.GetAlpha(ijk);

                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    Tensor3 dir = q[ijk] * stencil.Ct3(d);
                    size_t index = Index(ijk, d);
                    double cDotF = dir[1] * Fx[ijk] + dir[2] * Fy[ijk] + dir[3] * Fz[ijk];

                    double Gamma = eta[ijk] + kappa0[ijk] * E[ijk] + 3.0 * kappa1[ijk] * cDotF - I[index] * (kappaA[ijk] + kappa0[ijk]);
                    I[index] = I[index] + alpha * grid.dt * Gamma;
                }
            }
}
void Radiation::CollideStaticFluidBackwardEuler()
{
    PROFILE_FUNCTION();

    // int highestIteration = 0;
    // std::mutex mutex;
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                if (metric.InsideBH(grid.xyz(i, j, k)))
                    continue;
                size_t ijk = grid.Index(i, j, k);

                double alpha = metric.GetAlpha(ijk);
                double guessGammaLinear = 1.0 + alpha * grid.dt * (kappaA[ijk] + kappa0[ijk]);

                int cycles = 0;
                double diff = 1;
                while (abs(diff) > LAMBDA_ITTERATION_TOLERENCE && (cycles < MAX_LAMBDA_ITERATIONS))
                {
                    cycles++;
                    double guessE = E[ijk];
                    double guessFx = Fx[ijk];
                    double guessFy = Fy[ijk];
                    double guessFz = Fz[ijk];
                    double partOfGammaNoneLinear = eta[ijk] + kappa0[ijk] * guessE;

                    E[ijk] = 0.0;
                    Fx[ijk] = 0.0;
                    Fy[ijk] = 0.0;
                    Fz[ijk] = 0.0;
                    for (size_t d = 0; d < stencil.nDir; d++)
                    {
                        Tensor3 dir = q[ijk] * stencil.Ct3(d);
                        size_t index = Index(ijk, d);
                        double cDotGuessF = dir[1] * guessFx + dir[2] * guessFy + dir[3] * guessFz;

                        double guessGammaNoneLinear = partOfGammaNoneLinear + 3.0 * kappa1[ijk] * cDotGuessF;
                        Inew[index] = (I[index] + alpha * grid.dt * guessGammaNoneLinear) / guessGammaLinear;

                        double c = stencil.W(d) * Inew[index];
                        E[ijk] += c;
                        Fx[ijk] += c * dir[1];
                        Fy[ijk] += c * dir[2];
                        Fz[ijk] += c * dir[3];
                    }
                    double diffE = IntegerPow<2>((guessE - E[ijk]) / E[ijk]);
                    double diffFx = IntegerPow<2>((guessFx - Fx[ijk]) / Fx[ijk]);
                    double diffFy = IntegerPow<2>((guessFy - Fy[ijk]) / Fy[ijk]);
                    double diffFz = IntegerPow<2>((guessFz - Fz[ijk]) / Fz[ijk]);
                    diff = sqrt(diffE + diffFx + diffFy + diffFz);
                }
                //{
                //    std::lock_guard<std::mutex> lock(mutex);
                //    highestIteration = std::max(highestIteration, cycles);
                //}
            }
    // std::cout << " highestIteration = " << highestIteration << std::endl;
    std::swap(I, Inew);
}
void Radiation::CollideForwardEuler()
{
    PROFILE_FUNCTION();

    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                if (metric.InsideBH(grid.xyz(i, j, k)))
                    continue;
                size_t ijk = grid.Index(i, j, k);

                double alpha = metric.GetAlpha(ijk);
                Tensor3 u(0.0);                                                                          // fluid 3 velocity as seen by Eulerian observer
                double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ijk)));                          // Lorentz factor
                double uDotF = u[1] * Fx[ijk] + u[2] * Fy[ijk] + u[3] * Fz[ijk];                         // u_i F^i
                double uuDotP = u[1] * u[1] * Pxx[ijk] + u[2] * u[2] * Pyy[ijk] + u[3] * u[3] * Pzz[ijk] // u_i u_j P^ijk
                                + 2.0 * (u[1] * u[2] * Pxy[ijk] + u[1] * u[3] * Pxz[ijk] + u[2] * u[3] * Pyz[ijk]);
                double uDotPx = u[1] * Pxx[ijk] + u[2] * Pxy[ijk] + u[3] * Pxz[ijk]; // u_j P^ijk, i=1
                double uDotPy = u[1] * Pxy[ijk] + u[2] * Pyy[ijk] + u[3] * Pyz[ijk]; // u_j P^ijk, i=2
                double uDotPz = u[1] * Pxz[ijk] + u[2] * Pyz[ijk] + u[3] * Pzz[ijk]; // u_j P^ijk, i=3
                double fluidE = W * W * (E[ijk] - 2.0 * uDotF + uuDotP);
                double fluidFx = W * W * W * (2.0 * uDotF - E[ijk] - uuDotP) * u[1] + W * (Fx[ijk] - uDotPx);
                double fluidFy = W * W * W * (2.0 * uDotF - E[ijk] - uuDotP) * u[2] + W * (Fy[ijk] - uDotPy);
                double fluidFz = W * W * W * (2.0 * uDotF - E[ijk] - uuDotP) * u[3] + W * (Fz[ijk] - uDotPz);

                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    Tensor3 dir = q[ijk] * stencil.Ct3(d);
                    size_t index = Index(ijk, d);
                    double A = W * (1.0 - Tensor3::Dot(dir, u));
                    double cDotFluidF = dir[1] * fluidFx + dir[2] * fluidFy + dir[3] * fluidFz;

                    double Gamma = (eta[ijk] + kappa0[ijk] * fluidE + 3.0 * kappa1[ijk] * cDotFluidF) / (A * A * A) - A * I[index] * (kappaA[ijk] + kappa0[ijk]);
                    I[index] = std::max(I[index] + alpha * grid.dt * Gamma, 0.0);
                }
            }
}
void Radiation::CollideBackwardEuler()
{
    PROFILE_FUNCTION();

    int highestIteration = 0;
    // std::mutex mutex;
    PARALLEL_FOR(3)
    for (size_t k = HALO; k < grid.nz - HALO; k++)
        for (size_t j = HALO; j < grid.ny - HALO; j++)
            for (size_t i = HALO; i < grid.nx - HALO; i++)
            {
                if (metric.InsideBH(grid.xyz(i, j, k)))
                    continue;
                size_t ijk = grid.Index(i, j, k);

                double alpha = metric.GetAlpha(ijk);
                Tensor3 u(0.0);                                                 // fluid 2 velocity as seen by Eulerian observer
                double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ijk))); // Lorentz factor

                int cycles = 0;
                double diff = 1;
                while (abs(diff) > LAMBDA_ITTERATION_TOLERENCE && (cycles < MAX_LAMBDA_ITERATIONS))
                {
                    cycles++;
                    double guessE = E[ijk];
                    double guessFx = Fx[ijk];
                    double guessFy = Fy[ijk];
                    double guessFz = Fz[ijk];
                    double guessPxx = Pxx[ijk];
                    double guessPxy = Pxy[ijk];
                    double guessPxz = Pxz[ijk];
                    double guessPyy = Pyy[ijk];
                    double guessPyz = Pyz[ijk];
                    double guessPzz = Pzz[ijk];

                    double uDotF = u[1] * guessFx + u[2] * guessFy + u[3] * guessFz;                         // u_i F^i
                    double uuDotP = u[1] * u[1] * guessPxx + u[2] * u[2] * guessPyy + u[3] * u[3] * guessPzz // u_i u_j P^ijk
                                    + 2.0 * (u[1] * u[2] * guessPxy + u[1] * u[3] * guessPxz + u[2] * u[3] * guessPyz);
                    double uDotPx = u[1] * guessPxx + u[2] * guessPxy + u[3] * guessPxz; // u_j P^ijk, i=1
                    double uDotPy = u[1] * guessPxy + u[2] * guessPyy + u[3] * guessPyz; // u_j P^ijk, i=2
                    double uDotPz = u[1] * guessPxz + u[2] * guessPyz + u[3] * guessPzz; // u_j P^ijk, i=3
                    double guessFluidE = W * W * (guessE - 2.0 * uDotF + uuDotP);
                    double guessFluidFx = W * W * W * (2.0 * uDotF - guessE - uuDotP) * u[1] + W * (guessFx - uDotPx);
                    double guessFluidFy = W * W * W * (2.0 * uDotF - guessE - uuDotP) * u[2] + W * (guessFy - uDotPy);
                    double guessFluidFz = W * W * W * (2.0 * uDotF - guessE - uuDotP) * u[3] + W * (guessFz - uDotPz);

                    E[ijk] = 0.0;
                    Fx[ijk] = 0.0;
                    Fy[ijk] = 0.0;
                    Fz[ijk] = 0.0;
                    Pxx[ijk] = 0.0;
                    Pxy[ijk] = 0.0;
                    Pxz[ijk] = 0.0;
                    Pyy[ijk] = 0.0;
                    Pyz[ijk] = 0.0;
                    Pzz[ijk] = 0.0;
                    for (size_t d = 0; d < stencil.nDir; d++)
                    {
                        Tensor3 dir = q[ijk] * stencil.Ct3(d);
                        size_t index = Index(ijk, d);
                        double A = W * (1.0 - Tensor3::Dot(dir, u));
                        double cDotGuessFluidF = dir[1] * guessFluidFx + dir[2] * guessFluidFy + dir[3] * guessFluidFz;

                        double GuessGammaNoneLinear = (eta[ijk] + kappa0[ijk] * guessFluidE + 3.0 * kappa1[ijk] * cDotGuessFluidF) / (A * A * A);
                        double GuessGammaLinear = 1.0 + alpha * grid.dt * A * (kappaA[ijk] + kappa0[ijk]);
                        Inew[index] = (I[index] + alpha * grid.dt * GuessGammaNoneLinear) / GuessGammaLinear;

                        double c = stencil.W(d) * Inew[index];
                        E[ijk] += c;
                        Fx[ijk] += c * dir[1];
                        Fy[ijk] += c * dir[2];
                        Fz[ijk] += c * dir[3];
                        Pxx[ijk] += c * dir[1] * dir[1];
                        Pxy[ijk] += c * dir[1] * dir[2];
                        Pxz[ijk] += c * dir[1] * dir[3];
                        Pyy[ijk] += c * dir[2] * dir[2];
                        Pyz[ijk] += c * dir[2] * dir[3];
                        Pzz[ijk] += c * dir[3] * dir[3];
                    }
                    double diffE = abs((guessE - E[ijk]) / E[ijk]);
                    double diffFx = abs((guessFx - Fx[ijk]) / Fx[ijk]);
                    double diffFy = abs((guessFy - Fy[ijk]) / Fy[ijk]);
                    double diffFz = abs((guessFz - Fz[ijk]) / Fz[ijk]);
                    double diffPxx = abs((guessPxx - Pxx[ijk]) / Pxx[ijk]);
                    double diffPxy = abs((guessPxy - Pxy[ijk]) / Pxy[ijk]);
                    double diffPxz = abs((guessPxz - Pxz[ijk]) / Pxz[ijk]);
                    double diffPyy = abs((guessPyy - Pyy[ijk]) / Pyy[ijk]);
                    double diffPyz = abs((guessPyz - Pyz[ijk]) / Pyz[ijk]);
                    double diffPzz = abs((guessPzz - Pzz[ijk]) / Pzz[ijk]);
                    diff = diffE + diffFx + diffFy + diffFz + diffPxx + diffPxy + diffPxz + diffPyy + diffPyz + diffPzz;
                }
                //{
                //    std::lock_guard<std::mutex> lock(mutex);
                //    highestIteration = std::max(highestIteration, cycles);
                //}
            }
    // std::cout << " highestIteration = " << highestIteration << std::endl;
    std::swap(I, Inew);
}

void Radiation::TakePicture()
{
    // Note:
    // ijk is camera index
    // i,j,k are grid indexes and not related to ijk.
    PROFILE_FUNCTION();
    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < camera.pixelCount; ijk++)
    {
        Coord pixel = camera.xyz(ijk); // xyz coord of pixel in world space.
        if (grid.OutsideDomain(pixel))
        {
            camera.image[ijk] = 0;
            continue;
        }

        // We want to measure only the intensities that move orthogonal through the camera plane,
        // meaning the light velocity is the opposite of the camera normal vector.
        Tensor3 lookDir = camera.lookDirection;
        Tensor4 uLF(1, -lookDir[1], -lookDir[2], -lookDir[3]);
        uLF = NullNormalize(uLF, metric.GetMetric_ll(pixel));
        Tensor3 vIF = Vec3ObservedByEulObs<LF, IF>(uLF, pixel, metric);

        // Get 8 nearest Grid Points:
        double iTemp = grid.i(pixel[1]);
        double jTemp = grid.j(pixel[2]);
        double kTemp = grid.k(pixel[3]);
        size_t i0 = std::floor(iTemp);
        size_t i1 = i0 + 1;
        size_t j0 = std::floor(jTemp);
        size_t j1 = j0 + 1;
        size_t k0 = std::floor(kTemp);
        size_t k1 = k0 + 1;

        // Intensity interpolation:
        // double alpha = metric.GetAlpha(pixel); // removed to get intensity as seen by observer infinitly far away.
        double intensityAt_i0j0k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0, j0, k0))) * IntensityAt(grid.Index(i0, j0, k0), vIF);
        double intensityAt_i0j0k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0, j0, k1))) * IntensityAt(grid.Index(i0, j0, k1), vIF);
        double intensityAt_i0j1k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0, j1, k0))) * IntensityAt(grid.Index(i0, j1, k0), vIF);
        double intensityAt_i0j1k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0, j1, k1))) * IntensityAt(grid.Index(i0, j1, k1), vIF);
        double intensityAt_i1j0k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1, j0, k0))) * IntensityAt(grid.Index(i1, j0, k0), vIF);
        double intensityAt_i1j0k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1, j0, k1))) * IntensityAt(grid.Index(i1, j0, k1), vIF);
        double intensityAt_i1j1k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1, j1, k0))) * IntensityAt(grid.Index(i1, j1, k0), vIF);
        double intensityAt_i1j1k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1, j1, k1))) * IntensityAt(grid.Index(i1, j1, k1), vIF);

        camera.image[ijk] = TrilinearInterpolation(iTemp - i0, jTemp - j0, kTemp - k0,
                                                   intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
                                                   intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
    }
}

void Radiation::RunSimulation()
{
    // Initialize Profiler:
    Profiler::Session &session = Profiler::Session::Get();
    session.Start(config.name, "output/" + config.name + "/profileResults.json");

    // -------------------- Initialization --------------------
    LoadInitialData();
    UpdateSphericalHarmonicsCoefficients();

    int timeSteps = ceil(config.simTime / grid.dt);
    config.simTime = timeSteps * grid.dt;
    logger.SetValues(config.name, config.simTime);

    // Initial data output:
    if (config.printToTerminal)
    {
        std::cout << " nx           = " << grid.nx << "\n";
        std::cout << " ny           = " << grid.ny << "\n";
        std::cout << " nz           = " << grid.nz << "\n";
        std::cout << " nDir         = " << stencil.nDir << "\n";
        std::cout << " nSphHarm     = " << streamingStencil.nDir << "\n";
        std::cout << " simTime      = " << config.simTime << "\n";
        std::cout << " writePeriod  = " << config.writePeriod << "\n";
        std::cout << " dx           = " << grid.dx << "\n";
        std::cout << " dy           = " << grid.dy << "\n";
        std::cout << " dz           = " << grid.dz << "\n";
        std::cout << " dt           = " << grid.dt << "\n";
        std::cout << " timeSteps    = " << logger.timeSteps << "\n";
        std::cout << " filesToWrite = " << std::floor(config.simTime / config.writePeriod) << std::endl;
    }
    // --------------------------------------------------------

    // ----------------- Main simulation Loop -----------------
    {
        PROFILE_SCOPE("Total Time");
        double currentTime = config.t0;
        double timeSinceLastFrame = 0;

        // Save initial data:
        if (config.writeData)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, Fz_LF, logger.directoryPath + "/Moments/");
            timeSinceLastFrame = 0;
        }

        for (int n = 0; n < logger.timeSteps; n++)
        {
            if (config.printToTerminal)
                std::cout << "\n"
                          << n << "," << Format(currentTime, 4) << "," << std::flush;

            // Stream:
            switch (config.streamingType)
            {
            case (StreamingType::FlatFixed):
                StreamFlatFixed();
                break;
            case (StreamingType::FlatAdaptive):
                UpdateQuaternions();
                StreamFlatAdaptive();
                break;
            case (StreamingType::CurvedFixed):
                StreamCurvedFixed();
                break;
            case (StreamingType::CurvedAdaptive):
                UpdateQuaternions();
                StreamCurvedAdaptive();
                break;
            }

            // Collide:
            ComputeMomentsIF();
            // CollideStaticFluidForwardEuler();
            CollideStaticFluidBackwardEuler();
            // CollideForwardEuler();
            // CollideBackwardEuler();

            currentTime += grid.dt;
            timeSinceLastFrame += grid.dt;

            // Save data:
            if (config.writeData && timeSinceLastFrame >= config.writePeriod)
            {
                ComputeMomentsLF();
                grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, Fz_LF, logger.directoryPath + "/Moments/");
                timeSinceLastFrame = 0;
            }

            // Update other stuff:
            if (config.updateSphericalHarmonics)
                UpdateSphericalHarmonicsCoefficients();
            if (config.keepSourceNodesActive)
                LoadInitialData();
        }

        // Save final data:
        if (config.writeData)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, Fz_LF, logger.directoryPath + "/Moments/");
            timeSinceLastFrame = 0;
        }
    }
    // --------------------------------------------------------
    // Terminate profiler:
    session.End();

    // ---------------------- Termination ---------------------
    std::vector<std::string> names = session.GetAllFunctionNames();
    if (config.printToTerminal)
        std::cout << std::endl;
    for (int i = 0; i < names.size(); i++)
    {
        if (config.printToTerminal)
            session.PrintFunctionDuration(names[i]);
        logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
    }
    // --------------------------------------------------------
}