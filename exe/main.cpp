#include <iostream>
#include "../src/Radiation.h"
using namespace std;

// Macros:
#define WRITE_DATA true
#define PRINT_SETUP true
#define PRINT_PROGRESS true
#define PRINT_RESULTS true
#define SAVE_ID false

Logger SphereWave(LebedevStencil stencil, StreamingType streamingType, double cfl)
{
    // Grid, Metric, Stencil:
    size_t nx = 181;
    size_t ny = 181;
    size_t nz = 181;
    Coord start(-0.9, -0.9, -0.9);
    Coord end(0.9, 0.9, 0.9);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Config:
    Config config =
        {
            .name = "Sphere Wave 3d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 0.7,
            .writePeriod = 1.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    // Initial Data:
    double emissionRadius = 0.05;
    for(size_t k = 0; k < grid.nz; k++)
        for(size_t j = 0; j < grid.ny; j++)
            for(size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i,j,k);
                Coord xyz = grid.xyz(i,j,k);
                double r = xyz.EuklNorm();
                // Fluid:
                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                radiation.ux[ijk] = 0;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;
                if (r < emissionRadius)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    // Initial data given by moments:
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                    // Initial data given by intensities:
                    for(int d=0; d<stencil.nDir; d++)
                        radiation.initialI[radiation.Index(ijk, d)] = 1;
                }
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void SphereWaveAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) SphereWave(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 1) SphereWave(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) SphereWave(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 3) SphereWave(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 4) SphereWave(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 5) SphereWave(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 6) SphereWave(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 7) SphereWave(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 8) SphereWave(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);
}

Logger Shadow(LebedevStencil stencil, StreamingType streamingType, double cfl)
{
    // Grid, Metric, Stencil:
    size_t nx = 191;
    size_t ny = 191;
    size_t nz = 191;
    Coord start(-0.2, -0.2, -0.2);
    Coord end(1.7, 1.7, 1.7);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Config:
    Config config =
        {
            .name = "Shadow 3d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 1.5,
            .writePeriod = 2.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    // Initial Data:
    Coord planetPos(0.75, 0.75, 0.0);
    double sunRadius = 0.10;
    double planetRadius = 0.25;
    for(size_t k = 0; k < grid.nz; k++)
        for(size_t j = 0; j < grid.ny; j++)
            for(size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i,j,k);
                Coord xyz = grid.xyz(i,j,k);
                double r = xyz.EuklNorm();
                // Fluid:
                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                radiation.ux[ijk] = 0;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;
                if (r < sunRadius)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
                double dist = (planetPos - xyz).EuklNorm();
                if (dist <= planetRadius)
                    radiation.kappaA[ijk] = 1e10;
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void ShadowAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) Shadow(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 1) Shadow(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) Shadow(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 3) Shadow(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 4) Shadow(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 5) Shadow(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 6) Shadow(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 7) Shadow(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 8) Shadow(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);
}

Logger Star(LebedevStencil stencil, StreamingType streamingType, double cfl, double kappaA)
{
    // Grid, Metric, Stencil:
    size_t nx = 161;
    size_t ny = 161;
    size_t nz = 161;
    Coord start(-4, -4, -4);
    Coord end(4, 4, 4);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Config:
    Config config =
        {
            .name = "Star 3d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl" + Format(kappaA, 0) + "kappaA",
            .t0 = 0,
            .simTime = 10,
            .writePeriod = 11,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    // Initial Data (E0 and F1 assume starRadius = 1):
    double starRadius = 1;
    double E0 = 1.0 - exp(-kappaA);
    double F1 = (kappaA < 100) ? ((2.0 * kappaA * kappaA - 1.0) * exp(2.0 * kappaA) + 2.0 * kappaA + 1.0) / (8.0 * kappaA * kappaA * exp(2.0 * kappaA)) : 0.25;
    for(size_t k = 0; k < grid.nz; k++)
        for(size_t j = 0; j < grid.ny; j++)
            for(size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i,j,k);
                Coord xyz = grid.xyz(i,j,k);
                double r = xyz.EuklNorm();
                // Fluid:
                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.ux[ijk] = 0;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;
                if (r <= starRadius + grid.dx / 2.0)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.kappaA[ijk] = radiation.eta[ijk] = kappaA / starRadius;
                    radiation.initialE_LF[ijk] = E0;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
                else
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.kappaA[ijk] = radiation.eta[ijk] = 0;
                    radiation.initialE_LF[ijk]  = E0 / (r * r);
                    radiation.initialFx_LF[ijk] = F1 / (r * r) * xyz[1] / r;
                    radiation.initialFy_LF[ijk] = F1 / (r * r) * xyz[2] / r;
                    radiation.initialFz_LF[ijk] = F1 / (r * r) * xyz[3] / r;
                }
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void StarAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) Star(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl, 1);
    if(n == 1) Star(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl, 10);
    if(n == 2) Star(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl, 1e10);
    
    if(n == 3) Star(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl, 1);
    if(n == 4) Star(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl, 10);
    if(n == 5) Star(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl, 1e10);
    
    if(n == 6) Star(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl, 1);
    if(n == 7) Star(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl, 10);
    if(n == 8) Star(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl, 1e10);
}

Logger BeamCrossing(LebedevStencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 201;
    size_t ny = 101;
    size_t nz = 50;
    Coord start(-0.5, -0.25, -0.125);
    Coord end(0.5, 0.25, 0.125);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;
    
    InitialDataType initialDataType = InitialDataType::Moments;
    if (streamingType == StreamingType::FlatFixed)
        initialDataType = InitialDataType::Intensities;

    // Config:
    Config config =
        {
            .name = "Beam Crossing 3d/" + StreamingName(streamingType) + "_" + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 0.85,
            .writePeriod = 1,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = initialDataType,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    // Initial Data:
    Tensor3 dir0 = Tensor3(0.2,  0.1, 0.0).EuklNormalized();
    Tensor3 dir1 = Tensor3(0.2, -0.1, 0.0).EuklNormalized();
    
    // Find nearest directions in stencil:
    float dist0 = 1e10;
    float dist1 = 1e10;
    float d0 = -1;
    float d1 = -1;
    for(int d=0; d<stencil.nDir; d++)
    {
        float dist = (dir0 - stencil.Ct3(d)).EuklNorm();
        if(dist < dist0)
        {
            dist0 = dist;
            d0 = d;
        }
        dist = (dir1 - stencil.Ct3(d)).EuklNorm();
        if(dist < dist1)
        {
            dist1 = dist;
            d1 = d;
        }
    }

    // Initial Data:
    PARALLEL_FOR(3)
    for (size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];
                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;

                // Beam 0, from bottom to top:
                if ( -0.45 - grid.dx < x && x <= -0.45 && -0.20 < y && y < -0.15)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = dir0[1];
                    radiation.initialFy_LF[ijk] = dir0[2];
                    // Only used for FlatFixed streaming:
                    radiation.initialI[radiation.Index(ijk, d0)] = 1.0 / stencil.W(d0);
                }
                // Beam 1, from bottom to top:
                if ( -0.45 - grid.dx < x && x <= -0.45 && 0.15 < y && y < 0.20)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = dir1[1];
                    radiation.initialFy_LF[ijk] = dir1[2];
                    // Only used for FlatFixed streaming:
                    radiation.initialI[radiation.Index(ijk, d1)] = 1.0 / stencil.W(d1);
                }
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void BeamCrossingAnalysis(int n)
{
    // All missing
    double cfl = 0.9;
    // if(n == 0) BeamCrossing(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    
    if(n == 1) BeamCrossing(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) BeamCrossing(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 3) BeamCrossing(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 4) BeamCrossing(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 5) BeamCrossing(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 6) BeamCrossing(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);
}

Logger Diffusion(LebedevStencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl, double correctionFactor)
{
    // Create Radiation object:
    size_t nx = 101;
    size_t ny = 101;
    size_t nz = 101;
    Coord start(-0.5, -0.5, -0.5);
    Coord end(0.5, 0.5, 0.5);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (3.0 * kappa0) * (1.0 + correctionFactor * PE);

    // Config:
    double t0 = 1;
    std::string name = "Diffusion 3d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string(nz) + "nz " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE " + Format(correctionFactor, 3);
    Config config =
        {
            .name = name,
            .t0 = t0,
            .simTime = 2.0,
            .writePeriod = 1.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    for (size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];
                double r = xyz.EuklNorm();
                radiation.kappa0[ijk] = kappa0;
                radiation.kappa1[ijk] = kappa1;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                radiation.isInitialGridPoint[ijk] = true;
                radiation.ux[ijk] = 0.0;
                radiation.uy[ijk] = 0.0;
                radiation.uz[ijk] = 0.0;

                double E = pow(1.0 / t0, 3.0 / 2.0) * exp(-r * r / (4.0 * D * t0));
                radiation.initialE_LF[ijk] = E;
                radiation.initialFx_LF[ijk] = (x * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
                radiation.initialFy_LF[ijk] = (y * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
                radiation.initialFz_LF[ijk] = (z * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void DiffusionAnalysis(int n)
{
    double lambda = 0;
    // double cfl = 0.5;
    double cfl = 0.9;
    double correctionFactor = 0.75;
    // Refined stencil + adaptive streaming for kappa=10^2 and kappa=10^5:
    if (n == 0) Diffusion(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   ,    100.0, lambda, cfl, correctionFactor);
    if (n == 1) Diffusion(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , 100000.0, lambda, cfl, correctionFactor);
    if (n == 2) Diffusion(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive,    100.0, lambda, cfl, correctionFactor);
    if (n == 3) Diffusion(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, 100000.0, lambda, cfl, correctionFactor);
}

Logger MovingDiffusion(LebedevStencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl, double correctionFactor, double ux)
{
    // Create Radiation object:
    size_t nx = 2 * 150 + 1;
    size_t ny = 2 *  50 + 1;
    size_t nz = 2 *  50 + 1;
    Coord start(-1.0, -0.5, -0.5);
    Coord end(2.0, 0.5, 0.5);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (3.0 * kappa0) * (1.0 + correctionFactor * PE);

    // Account for time dilation:
    // We use t=1,2,3 when looking at the system from the lab frame (fluid is moving light gets dragged by the fluid)
    // We use t=(1,2,3)/gamma(ux) when looking at the system from the fluid frame (fluid is at rest)
    double simTime = 2.0;
    if (ux == 0)
    {
        double gamma = 1.0 / sqrt(1.0 - 0.5 * 0.5);
        simTime = simTime / gamma;
    }

    // Config:
    double t0 = 1;
    std::string name = "Moving Diffusion 3d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string(nz) + "nz " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE " + Format(ux, 3) + "ux";
    Config config =
        {
            .name = name,
            .t0 = t0,
            .simTime = simTime,
            .writePeriod = simTime / 2.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    for (size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double r = xyz.EuklNorm();
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];
                radiation.kappa0[ijk] = kappa0;
                radiation.kappa1[ijk] = kappa1;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                radiation.isInitialGridPoint[ijk] = true;
                radiation.ux[ijk] = ux;
                radiation.uy[ijk] = 0.0;
                radiation.uz[ijk] = 0.0;

                double E = pow(1.0 / t0, 3.0 / 2.0) * exp(-r * r / (4.0 * D * t0));
                double Fx = (x * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
                double Fy = (y * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
                double Fz = (z * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));

                // Lab Frame Moments:
                Tensor4x4 boost = BoostMatrix(Tensor3(-ux, 0, 0));
                Tensor4x4 Tff = Tensor4x4(E, Fx, Fy, Fz, Fx, 0, 0, 0, Fy, 0, 0, 0, Fz, 0, 0, 0);
                Tensor4x4 Tlf(0);
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        for (int I = 0; I < 4; I++)
                            for (int J = 0; J < 4; J++)
                                Tlf[{i, j}] += Tff[{I, J}] * boost[{i, I}] * boost[{j, J}];
                            
                // Use FF E for both so that the scaling is the same
                radiation.initialE_LF[ijk] = E;//Tlf[{0, 0}];
                radiation.initialFx_LF[ijk] = Tlf[{1, 0}];
                radiation.initialFy_LF[ijk] = Tlf[{2, 0}];
                radiation.initialFz_LF[ijk] = Tlf[{3, 0}];
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void MovingDiffusionAnalysis(int n)
{
    double lambda = 0;
    double cfl = 0.5;
    double correctionFactor = 0.75;
    // Note:
    // This test only works with stencils that support a minimum fluxMax value
    // This minimum value is yet to be determined.

    if (n == 4) MovingDiffusion(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , 1000.0, lambda, cfl, correctionFactor, 0.0);
    if (n == 5) MovingDiffusion(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , 1000.0, lambda, cfl, correctionFactor, 0.5);
    if (n == 6) MovingDiffusion(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, 1000.0, lambda, cfl, correctionFactor, 0.0);
    if (n == 7) MovingDiffusion(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, 1000.0, lambda, cfl, correctionFactor, 0.5);
}

Logger CurvedBeam(LebedevStencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 201;
    size_t ny = 161;
    size_t nz =  21;
    Coord start(0, 0, -0.25);
    Coord end(5, 4, 0.25);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Config:
    Config config =
        {
            .name = "Curved Beam 3d/" + metric.Name() + " " + stencil.name + " " + std::to_string(nx) + "nx" + std::to_string(ny) + "ny" + std::to_string(nz) + "nz" + Format(cfl, 2) + "cfl " + StreamingName(streamingType),
            .simTime = 10.0,
            .writePeriod = 11.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, interpGrid, camera, config);

    // Initial Data:
    PARALLEL_FOR(3)
    for (size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];
                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                if (x <= 1.1 * grid.dx && 3.00 <= y && y <= 3.50)
                {
                    Tensor4 uLF(1, 1, 0, 0);
                    uLF = NullNormalize(uLF, metric.GetMetric_ll(ijk));
                    Tensor3 vLF = Vec3ObservedByEulObs<LF, LF>(uLF, xyz, metric);

                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 10 * vLF[1];
                    radiation.initialFy_LF[ijk] = 10 * vLF[2];
                    radiation.initialFz_LF[ijk] = 10 * vLF[3];
                }
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void CurvedBeamAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) CurvedBeam(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl);
    if(n == 1) CurvedBeam(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl);
    if(n == 2) CurvedBeam(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::CurvedAdaptive, cfl);

    if(n == 3) CurvedBeam(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl);
    if(n == 4) CurvedBeam(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl);
    if(n == 5) CurvedBeam(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::CurvedAdaptive, cfl);

    if(n == 6) CurvedBeam(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl);
    if(n == 7) CurvedBeam(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl);
    if(n == 8) CurvedBeam(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::CurvedAdaptive, cfl);
}

Logger ThinHalfDisk(LebedevStencil stencil, StreamingType streamingType, double cfl, int resolutionScale)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    size_t nx = resolutionScale * 28 + 1;
    size_t ny = resolutionScale * 16 + 1;
    size_t nz = resolutionScale * 32 + 1;
    Coord start(-14, -1, -18);
    Coord   end( 14, 15,  14);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, m, a);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    size_t resX = 400;
    size_t resY = 400;
    size_t width = 28;
    size_t height = 28;
    Coord position(0,0,-16);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 10 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Config:
    Config config =
        {
            .name = "Thin Half Disk 3d/" + metric.Name() + " " + stencil.name + " " + std::to_string(nx) + "nx" + std::to_string(ny) + "ny" + std::to_string(nz) + "nz" + Format(cfl, 2) + "cfl " + StreamingName(streamingType),
            .simTime = 50.0,
            .writePeriod = 10.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = true,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

    // Initial Data:
    double kappaA = 1.0;
    double E0 = 1.0 - exp(-kappaA);
    #pragma omp parallel for
    for(size_t k=0; k<grid.nz; k++)
        for(size_t j=0; j<grid.ny; j++)
            for(size_t i=0; i<grid.nx; i++)
            {
                size_t ijk = grid.Index(i,j,k);
                Coord xyz = grid.xyz(i,j,k);
                double radius = xyz.EuklNorm();
                double phi = xyz.Phi();

                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                radiation.ux[ijk] = 0;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;

                //// Disk: with velocity and star like emissivity (use config.keepSourceNodesActive = false)
                //double v = sqrt(m / radius);
                //if(diskInner <= radius && radius <= diskOuter && abs(xyz[2]) < 0.9 * grid.dy)
                //{
                //    radiation.isInitialGridPoint[ijk] = true;
                //    radiation.kappaA[ijk] = radiation.eta[ijk] = kappaA;
                //    radiation.initialE_LF[ijk] = E0;
                //    radiation.initialFx_LF[ijk] = -v * MySin(phi);
                //    radiation.initialFy_LF[ijk] = 0;
                //    radiation.initialFz_LF[ijk] =  v * MyCos(phi);
                //}
                
                // Disk: simplistic (use config.keepSourceNodesActive = true)
                if(diskInner <= radius && radius <= diskOuter && abs(xyz[2]) < 0.9 * grid.dy)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void ThinHalfDiskAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) ThinHalfDisk(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 4);
    if(n == 1) ThinHalfDisk(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 4);
    if(n == 2) ThinHalfDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 4);
    if(n == 3) ThinHalfDisk(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 6);
    if(n == 4) ThinHalfDisk(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 6);
    if(n == 5) ThinHalfDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 6);
}

Logger ThinFullDisk(LebedevStencil stencil, StreamingType streamingType, double cfl, int resolutionScale)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    size_t nx = resolutionScale * 28 + 1;
    size_t ny = resolutionScale * 25 + 1;
    size_t nz = resolutionScale * 32 + 1;
    Coord start(-14, -10, -18);
    Coord   end( 14,  15,  14);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, m, a);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    size_t resX = 400;
    size_t resY = 600;
    size_t width = 28;
    size_t height = 42;
    Coord position(0,0,-16);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 10 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);


    // Config:
    Config config =
        {
            .name = "Thin Full Disk 3d/" + metric.Name() + " " + stencil.name + " " + std::to_string(nx) + "nx" + std::to_string(ny) + "ny" + std::to_string(nz) + "nz" + Format(cfl, 2) + "cfl " + StreamingName(streamingType),
            .simTime = 50.0,
            .writePeriod = 10.0,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = true,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

    // Initial Data:
    double kappaA = 1.0;
    double E0 = 1.0 - exp(-kappaA);
    #pragma omp parallel for
    for(size_t k=0; k<grid.nz; k++)
        for(size_t j=0; j<grid.ny; j++)
            for(size_t i=0; i<grid.nx; i++)
            {
                size_t ijk = grid.Index(i,j,k);
                Coord xyz = grid.xyz(i,j,k);
                double radius = xyz.EuklNorm();
                double phi = xyz.Phi();

                radiation.kappa0[ijk] = 0;
                radiation.kappa1[ijk] = 0;
                radiation.kappaA[ijk] = 0;
                radiation.eta[ijk] = 0;
                radiation.ux[ijk] = 0;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;

                //// Disk: with velocity and star like emissivity (use config.keepSourceNodesActive = false)
                //double v = sqrt(m / radius);
                //if(diskInner <= radius && radius <= diskOuter && abs(xyz[2]) < 0.9 * grid.dy)
                //{
                //    radiation.isInitialGridPoint[ijk] = true;
                //    radiation.kappaA[ijk] = radiation.eta[ijk] = kappaA;
                //    radiation.initialE_LF[ijk] = E0;
                //    radiation.initialFx_LF[ijk] = -v * MySin(phi);
                //    radiation.initialFy_LF[ijk] = 0;
                //    radiation.initialFz_LF[ijk] =  v * MyCos(phi);
                //}
                
                // Disk: simplistic (use config.keepSourceNodesActive = true)
                if(diskInner <= radius && radius <= diskOuter && abs(xyz[2]) < 0.9 * grid.dy)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
            }
    radiation.RunSimulation();
    return radiation.logger;
}
void ThinFullDiskAnalysis(int n)
{
    double cfl = 0.9;
    // if(n == 0) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 3);
    // if(n == 1) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 4);
    // if(n == 2) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 5);
    // if(n == 3) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 6);
    // if(n == 4) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 7);
    // if(n == 5) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 8);
    // if(n == 6) ThinFullDisk(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 9);
    
    if(n == 0) ThinFullDisk(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 3);
    if(n == 1) ThinFullDisk(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 4);
    if(n == 2) ThinFullDisk(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 5);
    if(n == 3) ThinFullDisk(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 6);
    if(n == 4) ThinFullDisk(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::CurvedFixed, cfl, 7);
}

void StencilAnalysis(int n)
{
    int order[] = {11, 15, 17, 19, 21, 23, 29, 31, 35, 41, 47};
    if (n == 0)
        for(int i=0; i<11; i++)
        {// All good
            LebedevStencil stencil = LebedevStencil(order[i], 0.00, 0.00, 0.00);
            stencil.Print();
        }
    if (n == 1)
        for(int i=0; i<11; i++)
        {// All good
            LebedevStencil stencil = LebedevStencil(order[i], 0.15, 0.00, 0.00);
            stencil.Print();
        }
    if (n == 2)
        for(int i=0; i<11; i++)
        {
            LebedevStencil stencil = LebedevStencil(order[i], 0.00, 0.15, 0.00);
            stencil.Print();
        }
}
void TestThinDiskSetup()
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Camera
    size_t resX = 40;
    size_t resY = 60;
    size_t width = 28;
    size_t height = 42;
    Coord position(0,0,-16);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 10 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Grid:
    size_t nx = 280;
    size_t ny = 250;
    size_t nz = 320;
    Coord start(-14, -10, -18);
    Coord   end( 14,  15,  14);
    Grid grid(nx, ny, nz, start, end);

    std::ofstream fileOut(OUTPUTDIR + "TestThinDiskSetup.csv");
    fileOut << "#nx=" << grid.nx << "\n";
    fileOut << "#ny=" << grid.ny << "\n";
    fileOut << "#nz=" << grid.nz << "\n";
    fileOut << "#startx=" << grid.startx << "\n";
    fileOut << "#starty=" << grid.starty << "\n";
    fileOut << "#startz=" << grid.startz << "\n";
    fileOut << "#endx=" << grid.endx << "\n";
    fileOut << "#endy=" << grid.endy << "\n";
    fileOut << "#endz=" << grid.endz << "\n";
    fileOut << "#x,y,z,r,g,b,a\n";

    // Accretion disk:
    for(size_t k=0; k<grid.nz; k++)
        for(size_t j=0; j<grid.ny; j++)
            for(size_t i=0; i<grid.nx; i++)
            {
                Coord xyz = grid.xyz(i,j,k);
                double radius = xyz.EuklNorm();
                
                if(diskInner <= radius && radius <= diskOuter && abs(xyz[2]) < 0.9 * grid.dy)
                {
                    fileOut << xyz[1] << "," << xyz[2] << "," << xyz[3] << ",0,0,0,0\n";
                }
            }

    // Camera plane:
    for(size_t j=0; j<resY; j++)
        for(size_t i=0; i<resX; i++)
        {
            size_t ij = camera.Index(i ,j);
            Coord x = camera.xyzWorld(i, j);
            fileOut << x[1] << "," << x[2] << "," << x[3] << ",1,1,1,1\n";
        }

    fileOut.close();
}

int main(int argc, char *argv[])
{
    int n = 1;
    if (argc > 1)
        n = atoi(argv[1]);

    // SphereWaveAnalysis(n);       // Done
    // ShadowAnalysis(n);           // Done
    // StarAnalysis(n);             // Done
    // BeamCrossingAnalysis(n);     // Done
    // DiffusionAnalysis(n);        // Done
    // MovingDiffusionAnalysis(n);  // Done
    // CurvedBeamAnalysis(n);       // Done
    // TestThinDiskSetup();         // Done
    // ThinHalfDiskAnalysis(n);     // Not needed for paper, just a quick test for Full Disk
    // ThinFullDiskAnalysis(n);     // Done
    // StencilAnalysis(n);          // Done
}