#include "Radiation.h"

Radiation::Radiation(Metric &metric, LebedevStencil &stencil, LebedevStencil &streamingStencil, InterpolationGrid &interpGrid, Camera &camera, Config config)
    : grid(metric.grid), metric(metric), stencil(stencil), streamingStencil(streamingStencil), interpGrid(interpGrid), camera(camera), config(config), logger(stencil, streamingStencil, metric)
{
    isInitialGridPoint = new bool[grid.nxyz]();
    initialI.resize(grid.nxyz * stencil.nDir);
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
    ux.resize(grid.nxyz);
    uy.resize(grid.nxyz);
    uz.resize(grid.nxyz);

    I.resize(grid.nxyz * stencil.nDir);
    Inew.resize(grid.nxyz * stencil.nDir);
    
    coefficientsS.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsX.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsY.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsZ.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsCx.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsCy.resize(grid.nxyz * streamingStencil.nCoefficients);
    coefficientsCz.resize(grid.nxyz * streamingStencil.nCoefficients);

    itterationCount.resize(grid.nxyz);

    // Initialize all quaternions to identity:
    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
    {
        q[ijk] = qNew[ijk] = glm::quat(1, 0, 0, 0);
        itterationCount[ijk] = 0;
    }

    // Initialize all spherical harmonic coefficinets to identity:
    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
    {
            double dataS[streamingStencil.nDir];
            double dataX[streamingStencil.nDir];
            double dataY[streamingStencil.nDir];
            double dataZ[streamingStencil.nDir];
            double dataCx[streamingStencil.nDir];
            double dataCy[streamingStencil.nDir];
            double dataCz[streamingStencil.nDir];
            for (size_t d = 0; d < streamingStencil.nDir; d++)
            {
                Coord xyz = grid.xyz(ijk);
                dataS[d] = 1.0;
                dataX[d] = xyz[1];
                dataY[d] = xyz[2];
                dataY[d] = xyz[3];
                dataCx[d] = stencil.Cx(d);
                dataCy[d] = stencil.Cy(d);
                dataCz[d] = stencil.Cz(d);
            }
            SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataS, &coefficientsS[HarmonicIndex(0, ijk)]);
            SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataX, &coefficientsX[HarmonicIndex(0, ijk)]);
            SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataY, &coefficientsY[HarmonicIndex(0, ijk)]);
            SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[HarmonicIndex(0, ijk)]);
            SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[HarmonicIndex(0, ijk)]);
    }
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
    
    if (config.initialDataType == InitialDataType::Intensities)
    {
        if (isAdaptiveStreaming)
            LoadInitialDataIntensitiesAdaptive();
        else
            LoadInitialDataIntensitiesFixed();
    }
    else if (config.initialDataType == InitialDataType::Moments)
    {
        if (isAdaptiveStreaming)
            LoadInitialDataMomentsAdaptive();
        else
            LoadInitialDataMomentsFixed();
    }
    else
        ExitOnError("Radiation.cpp: LoadInitialData() => Unknown initial data type.");
}
void Radiation::LoadInitialDataIntensitiesFixed()
{
    // When using a fixed stencil extra ghost directions have no real benefit.
    if (stencil.refinement0Threshold > 0.0 || stencil.refinement1Threshold > 0.0 || stencil.refinement2Threshold > 0.0)
        ExitOnError("Radiation.cpp: LoadInitialDataIntensitiesFixed() => Refined stencils do not support fixed streaming.");

    PARALLEL_FOR(1)
    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
    {
        if (!isInitialGridPoint[ijk])
            continue;

        for (size_t d = 0; d < stencil.nDir; d++)
            I[Index(ijk, d)] = initialI[Index(ijk, d)];
    }
}
void Radiation::LoadInitialDataMomentsFixed()
{
    // When using a fixed stencil extra ghost directions have no real benefit.
    if (stencil.refinement0Threshold > 0.0 || stencil.refinement1Threshold > 0.0 || stencil.refinement2Threshold > 0.0)
        ExitOnError("Radiation.cpp: LoadInitialDataMomentsFixed() => Refined stencils do not support fixed streaming.");

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

        // Relative flux magnitude, sigma, flux direction:
        Tensor3 initialFxyz_IF(initialFx_IF, initialFy_IF, initialFz_IF);
        double initialF_IF = initialFxyz_IF.EuklNorm();
        double relativeF_IF = initialF_IF / initialE_IF;
        double sigma = stencil.sigmaOfRelativeFlux.Evaluate(relativeF_IF);
        Tensor3 dirInitialF = (initialF_IF < MIN_FLUX_NORM) ? Tensor3(0, 0, 1) : Tensor3(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF, initialFz_IF / initialF_IF);

        for (size_t d = 0; d < stencil.nDir; d++)
            I[Index(ijk, d)] = Intensity(sigma, initialE_IF, dirInitialF, stencil.Ct3(d));
    }
}
void Radiation::LoadInitialDataIntensitiesAdaptive()
{
    ExitOnError("Radiation.cpp: LoadInitialDataIntensitiesAdaptive() => Adaptive stencils do not support initial data by intensities.");
}
void Radiation::LoadInitialDataMomentsAdaptive()
{
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

        // Relative flux magnitude, sigma, stencil quaternion (to):
        Tensor3 initialFxyz_IF(initialFx_IF, initialFy_IF, initialFz_IF);
        double initialF_IF = initialFxyz_IF.EuklNorm();
        double relativeF_IF = initialF_IF / initialE_IF;
        double sigma = stencil.sigmaOfRelativeFlux.Evaluate(relativeF_IF);
        glm::vec3 to = (initialF_IF < MIN_FLUX_NORM) ? glm::vec3(0, 0, 1) : glm::vec3(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF, initialFz_IF / initialF_IF);

        q[ijk] = glm::quat(from, to);
        for (size_t d = 0; d < stencil.nDir; d++)
            I[Index(ijk, d)] = Intensity(sigma, initialE_IF, stencil.Theta(d));
    }
}
// Weighting initial data with 1/(w*N) produces worse results.
//void Radiation::LoadInitialDataMomentsAdaptive()
//{
//    // Base stencil and intensities for initial data interpolation:
//    LebedevStencil stencilBase(stencil.nOrder);
//    RealBuffer IBase;
//    IBase.resize(stencilBase.nDir * grid.nxyz);
//    auto BaseIndex = [&](size_t ijk, size_t d)
//    { return d + ijk * stencilBase.nDir; };
//    
//    PARALLEL_FOR(1)
//    for (size_t ijk = 0; ijk < grid.nxyz; ijk++)
//    {
//        if (!isInitialGridPoint[ijk])
//            continue;
//
//        // Convert given LF initial data to IF:
//        Tensor4 initialDataIF = InitialDataLFtoIF(ijk);
//        double initialE_IF = initialDataIF[0];
//        double initialFx_IF = initialDataIF[1];
//        double initialFy_IF = initialDataIF[2];
//        double initialFz_IF = initialDataIF[3];
//
//        // Relative flux magnitude, sigma, stencil quaternion (to):
//        Tensor3 initialFxyz_IF(initialFx_IF, initialFy_IF, initialFz_IF);
//        double initialF_IF = initialFxyz_IF.EuklNorm();
//        double relativeF_IF = initialF_IF / initialE_IF;
//        double sigma = stencil.sigmaOfRelativeFlux.Evaluate(relativeF_IF);
//        glm::vec3 to = (initialF_IF < MIN_FLUX_NORM) ? glm::vec3(0, 0, 1) : glm::vec3(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF, initialFz_IF / initialF_IF);
//
//        // Calculated deweighted intensities on stencilBase:
//        q[ijk] = glm::quat(from, to);
//        for (size_t d = 0; d < stencilBase.nDir; d++)
//            IBase[BaseIndex(ijk, d)] = Intensity(sigma, initialE_IF, stencilBase.Theta(d)) / (stencilBase.W(d) * stencilBase.nDir);
//
//        // Interpolate to actual stencil:
//        for (size_t d = 0; d < stencil.nDir; d++)
//        {
//            Vector3 dir = stencil.Cv3(d);
//            std::tuple<std::vector<size_t>, std::vector<double>> neighboursAndWeights = stencilBase.VoronoiNeighboursAndWeights(dir);
//            std::span<const size_t> neighbours = std::get<0>(neighboursAndWeights);
//            std::span<const double> weights = std::get<1>(neighboursAndWeights);
//            double interpolatetValue = 0;
//            for (size_t p = 0; p < weights.size(); p++)
//                interpolatetValue += weights[p] * IBase[BaseIndex(ijk, neighbours[p])];
//            I[Index(ijk,d)] = interpolatetValue;
//        }
//    }
//}

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
    // i,j,k >= 1, thus i+a etc will never be negative.
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

void Radiation::Collide()
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
                double uu = ux[ijk] * ux[ijk] + uy[ijk] * uy[ijk] + uz[ijk] * uz[ijk];
                double lorentz = 1.0 / sqrt(1.0 - uu);

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

                    double uF = ux[ijk] * guessFx + uy[ijk] * guessFy + uz[ijk] * guessFz;    // u_i F^i
                    double uxP = ux[ijk] * guessPxx + uy[ijk] * guessPxy + uz[ijk] * guessPxz; // x component of u_i P^ij
                    double uyP = ux[ijk] * guessPxy + uy[ijk] * guessPyy + uz[ijk] * guessPyz; // y component of u_i P^ij
                    double uzP = ux[ijk] * guessPxz + uy[ijk] * guessPyz + uz[ijk] * guessPzz; // y component of u_i P^ij
                    double uuP = ux[ijk] * uxP + uy[ijk] * uyP + uz[ijk] * uzP; // u_i u_j P^ij
                    double guessFluidE = lorentz * lorentz * (guessE - 2.0 * uF + uuP);
                    double guessFluidFx = lorentz * lorentz * (guessFluidFx - lorentz * guessE * ux[ijk] - uxP + ux[ijk] * lorentz / (1.0 + lorentz) * ((2.0 * lorentz + 1.0) * uF - lorentz * uuP));
                    double guessFluidFy = lorentz * lorentz * (guessFluidFy - lorentz * guessE * uy[ijk] - uyP + uy[ijk] * lorentz / (1.0 + lorentz) * ((2.0 * lorentz + 1.0) * uF - lorentz * uuP));
                    double guessFluidFz = lorentz * lorentz * (guessFluidFz - lorentz * guessE * uz[ijk] - uzP + uz[ijk] * lorentz / (1.0 + lorentz) * ((2.0 * lorentz + 1.0) * uF - lorentz * uuP));

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
                        Tensor3 nIF = q[ijk] * stencil.Ct3(d);
                        size_t index = Index(ijk, d);
                        double un = ux[ijk] * nIF[1] + uy[ijk] * nIF[2] + uz[ijk] * nIF[3];
                        double A = lorentz * (1.0 - un);
                        double nF = nIF[1] * guessFx + nIF[2] * guessFy + nIF[3] * guessFz;

                        double M = kappa0[ijk] * guessFluidE + 3.0 * kappa1[ijk] * nF;
                        Inew[index] = (I[index] + alpha * grid.dt * (eta[ijk] + M) / (A * A * A)) / (1.0 + alpha * grid.dt * (kappaA[ijk] + kappa0[ijk]) * A);

                        double c = stencil.W(d) * Inew[index];
                        E[ijk] += c;
                        Fx[ijk] += c * nIF[1];
                        Fy[ijk] += c * nIF[2];
                        Fz[ijk] += c * nIF[3];
                        Pxx[ijk] += c * nIF[1] * nIF[1];
                        Pxy[ijk] += c * nIF[1] * nIF[2];
                        Pxz[ijk] += c * nIF[1] * nIF[3];
                        Pyy[ijk] += c * nIF[2] * nIF[2];
                        Pyz[ijk] += c * nIF[2] * nIF[3];
                        Pzz[ijk] += c * nIF[3] * nIF[3];
                    }

                    double diffE = IntegerPow<2>((guessE - E[ijk]) / E[ijk]);
                    double diffFx = IntegerPow<2>((guessFx - Fx[ijk]) / Fx[ijk]);
                    double diffFy = IntegerPow<2>((guessFy - Fy[ijk]) / Fy[ijk]);
                    double diffFz = IntegerPow<2>((guessFz - Fz[ijk]) / Fz[ijk]);
                    double diffPxx = IntegerPow<2>((guessPxx - Pxx[ijk]) / Pxx[ijk]);
                    double diffPxy = IntegerPow<2>((guessPxy - Pxy[ijk]) / Pxy[ijk]);
                    double diffPxz = IntegerPow<2>((guessPxz - Pxz[ijk]) / Pxz[ijk]);
                    double diffPyy = IntegerPow<2>((guessPyy - Pyy[ijk]) / Pyy[ijk]);
                    double diffPyz = IntegerPow<2>((guessPyz - Pyz[ijk]) / Pyz[ijk]);
                    double diffPzz = IntegerPow<2>((guessPzz - Pzz[ijk]) / Pzz[ijk]);
                    diff = sqrt(diffE + diffFx + diffFy + diffFz + diffPxx + diffPxy + diffPxz + diffPyy + diffPyz + diffPzz);
                }
                itterationCount[ijk] = cycles;
            }
    
    int maxItteration = 0;
    for(int ijk = 0; ijk < grid.nxyz; ijk++)
    {
        maxItteration = std::max(maxItteration, itterationCount[ijk]);
        averageItterationCount += itterationCount[ijk];
        itterationCount[ijk] = 0;
    }
    maxItterationCount = std::max(maxItterationCount, maxItteration);

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
        Coord pixelWorld = camera.xyzWorld(ijk); // xyz coord of pixel in world space.
        if (grid.OutsideDomain(pixelWorld))
        {
            camera.image[ijk] = 0;
            continue;
        }

        // We want to measure only the intensities that moves orthogonal through the camera plane,
        // meaning the light velocity is the opposite of the camera normal vector.
        Tensor3 lookDir = camera.lookDirection;
        Tensor4 uLF(1, -lookDir[1], -lookDir[2], -lookDir[3]);
        uLF = NullNormalize(uLF, metric.GetMetric_ll(pixelWorld));
        Tensor3 vIF = Vec3ObservedByEulObs<LF, IF>(uLF, pixelWorld, metric);

        // Get 8 nearest Grid Points:
        double iTemp = grid.i(pixelWorld[1]);
        double jTemp = grid.j(pixelWorld[2]);
        double kTemp = grid.k(pixelWorld[3]);
        size_t i0 = std::floor(iTemp);
        size_t i1 = i0 + 1;
        size_t j0 = std::floor(jTemp);
        size_t j1 = j0 + 1;
        size_t k0 = std::floor(kTemp);
        size_t k1 = k0 + 1;

        // Intensity interpolation:
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
    if (config.printSetup)
    {
        std::cout << " nx           = " << grid.nx << "\n";
        std::cout << " ny           = " << grid.ny << "\n";
        std::cout << " nz           = " << grid.nz << "\n";
        std::cout << " nDir         = " << stencil.nDir << "\n";
        std::cout << " sigmaMax     = " << stencil.sigmaMax << "\n";
        std::cout << " fluxMax      = " << stencil.relativeFluxMax << "\n";
        std::cout << " simTime      = " << config.simTime << "\n";
        std::cout << " nSphHarm     = " << streamingStencil.nDir << "\n";
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
        if (config.writeData && config.saveInitialData)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, Fz_LF, logger.directoryPath + "/Moments/");
            if (config.useCamera)
            {
                TakePicture();
                camera.WriteImagetoCsv(currentTime, logger.directoryPath + "/Images/");
            }
            timeSinceLastFrame = 0;
        }

        for (int n = 0; n < logger.timeSteps; n++)
        {
            if (config.printProgress)
                std::cout << "\n" << n << "," << Format(currentTime, 4) << "," << std::flush;

            // Update stuff:
            if (config.updateSphericalHarmonics)
                UpdateSphericalHarmonicsCoefficients();
            if (config.keepSourceNodesActive)
                LoadInitialData();

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
            Collide();

            currentTime += grid.dt;
            timeSinceLastFrame += grid.dt;

            // Save data:
            if (config.writeData && timeSinceLastFrame >= config.writePeriod - 1e-8)
            {
                ComputeMomentsIF();
                ComputeMomentsLF();
                grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, Fz_LF, logger.directoryPath + "/Moments/");
                if (config.useCamera)
                {
                    TakePicture();
                    camera.WriteImagetoCsv(currentTime, logger.directoryPath + "/Images/");
                }
                timeSinceLastFrame = 0;
            }
        }

        // Save final data:
        if (config.writeData && timeSinceLastFrame != 0)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, Fz_LF, logger.directoryPath + "/Moments/");
            if (config.useCamera)
            {
                TakePicture();
                camera.WriteImagetoCsv(currentTime, logger.directoryPath + "/Images/");
            }
            timeSinceLastFrame = 0;
        }
    }
    // --------------------------------------------------------
    // Terminate profiler:
    session.End();

    // ---------------------- Termination ---------------------
    if (config.printResults)
        std::cout << std::endl;
        
    logger.maxItterationCount = maxItterationCount;
    averageItterationCount /= (logger.timeSteps * (grid.nx - 2) * (grid.ny - 2) * (grid.nz - 2));
    logger.averageItterationCount = averageItterationCount;
    if (config.printResults)
    {
        std::cout << "Max Lambda Itteration Count     = " << maxItterationCount << std::endl;
        std::cout << "Average Lambda Itteration Count = " << averageItterationCount << std::endl;
    }

    std::vector<std::string> names = session.GetAllFunctionNames();
    double totalTime;
    double writingTime;
    for (int i = 0; i < names.size(); i++)
    {
        if (config.printResults)
            session.PrintFunctionDuration(names[i]);
        if (names[i] == "Total Time")
            totalTime = session.GetTotalTime(names[i]);
        if (names[i] == "void Grid::WriteFrametoCsv(float, const RealBuffer&, const RealBuffer&, const RealBuffer&, const RealBuffer&, std::string, std::string)")
            writingTime = session.GetTotalTime(names[i]);
        logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
    }
    if (config.printResults)
        std::cout << "Computation Time: " << totalTime - writingTime << "s" << std::endl;
    // --------------------------------------------------------
}