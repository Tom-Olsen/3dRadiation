#include <iostream>
#include "../src/Radiation.h"
using namespace std;

#define SAVE_ID false
#define PRINT_SETUP true
#define PRINT_PROGRESS true
#define PRINT_RESULTS true

void SphereWave(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Grid, Metric, Stencil:
    size_t nx = 180 + 1 + 2;
    size_t ny = 180 + 1 + 2;
    size_t nz = 180 + 1 + 2;
    double dx = 1.8 / (nx - 1.0 - 2.0);
    double dy = 1.8 / (ny - 1.0 - 2.0);
    double dz = 1.8 / (nz - 1.0 - 2.0);
    Coord start(-0.9 - dx, -0.9 - dy, -0.9 - dz);
    Coord end(0.9 + dx, 0.9 + dy, 0.9 + dz);
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
            .writeData = true,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Intensities,
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
                    // Initial data given by intensities:
                    for(int d=0; d<stencil.nDir; d++)
                        radiation.initialI[radiation.Index(ijk, d)] = 1;
                    glm::vec3 from = glm::vec3(0, 0, 1);
                    glm::vec3 to = (r > 1e-6) ? glm::vec3(xyz[1] / r, xyz[2] / r, xyz[3] / r) : from;
                    radiation.initialQ[ijk] = glm::quat(from, to);
                }
            }
    radiation.RunSimulation();
}
void SphereWaveAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) SphereWave(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 1) SphereWave(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) SphereWave(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 0) SphereWave(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 1) SphereWave(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) SphereWave(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 0) SphereWave(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    if(n == 1) SphereWave(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) SphereWave(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);
}

void Shadow(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Grid, Metric, Stencil:
    size_t nx = 190 + 1 + 2;
    size_t ny = 190 + 1 + 2;
    size_t nz = 190 + 1 + 2;
    double dx = 1.1 / (nx - 1.0 - 2.0);
    double dy = 1.1 / (ny - 1.0 - 2.0);
    double dz = 1.1 / (nz - 1.0 - 2.0);
    Coord start(-0.2 - dx, -0.2 - dy, -0.2 - dz);
    Coord end(1.7 + dx, 1.7 + dy, 1.7 + dz);
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
            .keepSourceNodesActive = false,
            .writeData = true,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Intensities,
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
                    radiation.kappa0[ijk] = radiation.eta[ijk] = 1;
                    radiation.isInitialGridPoint[ijk] = true;
                    //radiation.isInitialGridPoint[ijk] = true;
                    //// Initial data given by intensities:
                    //for(int d=0; d<stencil.nDir; d++)
                    //        radiation.initialI[radiation.Index(ijk, d)] = 1.0 / (stencil.W(d) * stencil.nDir);
                    //glm::vec3 from = glm::vec3(0, 0, 1);
                    //glm::vec3 to = (r > 1e-6) ? glm::vec3(xyz[1] / r, xyz[2] / r, xyz[3] / r) : from;
                    //radiation.initialQ[ijk] = glm::quat(from, to);
                }
            double dist = (planetPos - xyz).EuklNorm();
            if (dist <= planetRadius)
                radiation.kappaA[ijk] = 1e10;
            }
    radiation.RunSimulation();
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

void BeamCrossing(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 200 + 1 + 2;
    size_t ny = 100 + 1 + 2;
    size_t nz = 50 + 1 + 2;
    double dx = 4.0 / (nx - 1.0 - 2.0);
    double dy = 2.0 / (ny - 1.0 - 2.0);
    double dz = 1.0 / (nz - 1.0 - 2.0);
    Coord start(-0.5 - dx, -0.25 - dy, -0.125 - dz);
    Coord end(0.5 + dx, 0.25 + dy, 0.125 + dz);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
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
            .writeData = true,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = initialDataType,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

    // Initial Data:
    Tensor3 dir0 = Tensor3(0.3,  0.1, 0.0).EuklNormalized();
    Tensor3 dir1 = Tensor3(0.3, -0.1, 0.0).EuklNormalized();
    
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
                if ( -0.45 - dx < x && x <= -0.45 && -0.20 < y && y < -0.15)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = dir0[1];
                    radiation.initialFy_LF[ijk] = dir0[2];
                    // Only used for FlatFixed streaming:
                    radiation.initialI[radiation.Index(ijk, d0)] = 1.0 / stencil.W(d0);
                }
                // Beam 1, from bottom to top:
                if ( -0.45 - dx < x && x <= -0.45 && 0.15 < y && y < 0.20)
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
}
void BeamCrossingAnalysis(int n)
{
    // All missing
    double cfl = 0.9;
    if(n == 0) BeamCrossing(LebedevStencil(59, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , cfl);
    
    if(n == 1) BeamCrossing(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 2) BeamCrossing(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 3) BeamCrossing(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 4) BeamCrossing(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);

    if(n == 5) BeamCrossing(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive, cfl);
    if(n == 6) BeamCrossing(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, cfl);
    
}

void Diffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl, double correctionFactor)
{
    // Create Radiation object:
    size_t nx = 100 + 1 + 2;
    size_t ny = 100 + 1 + 2;
    size_t nz = 100 + 1 + 2;
    double dx = 1.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    double dz = 1.0 / (nz - 1.0 - 2.0);
    Coord start(-0.5 - dx, -0.5 - dy, -0.5 - dz);
    Coord end(0.5 + dx, 0.5 + dy, 0.5 + dz);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
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
            .writeData = true,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

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
}
void DiffusionAnalysis(int n)
{
    double lambda = 0;
    double cfl = 0.5;
    double correctionFactor = 0.75;
    // Refined stencil + adaptive streaming for kappa=10^2 and kappa=10^5:
    if (n == 0) Diffusion(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   ,    100.0, lambda, cfl, correctionFactor);
    if (n == 1) Diffusion(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::FlatFixed   , 100000.0, lambda, cfl, correctionFactor);
    if (n == 2) Diffusion(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive,    100.0, lambda, cfl, correctionFactor);
    if (n == 3) Diffusion(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::FlatAdaptive, 100000.0, lambda, cfl, correctionFactor);
}

void MovingDiffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl, double correctionFactor, double ux)
{
    // Create Radiation object:
    size_t nx = 2 * 150 + 1 + 2;
    size_t ny = 2 * 50 + 1 + 2;
    size_t nz = 2 * 50 + 1 + 2;
    double dx = 3.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    double dz = 1.0 / (nz - 1.0 - 2.0);
    Coord start(-1.0 - dx, -0.5 - dy, -0.5 - dz);
    Coord end(2.0 + dx, 0.5 + dy, 0.5 + dz);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
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
            .writeData = true,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

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

void CurvedBeam(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 200 + 1 + 2;
    size_t ny = 160 + 1 + 2;
    size_t nz = 20 + 1 + 2;
    double dx = 5.0 / (nx - 1.0 - 2.0);
    double dy = 4.0 / (ny - 1.0 - 2.0);
    double dz = 0.5 / (nz - 1.0 - 2.0);
    Coord start(0 - dx, 0 - dy, -0.25 - dz);
    Coord end(5 + dx, 4 + dy, 0.25 + dz);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
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
            .writeData = true,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .useCamera = false,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

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
}
void CurvedBeamAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) CurvedBeam(LebedevStencil(29, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl);
    if(n == 1) CurvedBeam(LebedevStencil(29, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl);
    if(n == 2) CurvedBeam(LebedevStencil(29, 0.00, 0.15, 0.00), StreamingType::CurvedAdaptive, cfl);

    if(n == 0) CurvedBeam(LebedevStencil(35, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl);
    if(n == 1) CurvedBeam(LebedevStencil(35, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl);
    if(n == 2) CurvedBeam(LebedevStencil(35, 0.00, 0.15, 0.00), StreamingType::CurvedAdaptive, cfl);

    if(n == 0) CurvedBeam(LebedevStencil(41, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl);
    if(n == 1) CurvedBeam(LebedevStencil(41, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl);
    if(n == 2) CurvedBeam(LebedevStencil(41, 0.00, 0.15, 0.00), StreamingType::CurvedAdaptive, cfl);
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

/*
void ThinDisk1(size_t nx, size_t ny, size_t nz, Stencil stencil, int sigma, int simTime, StreamingType streamingType, string comment)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    Coord start(-14,-14,-1);
    Coord   end( 14, 22, 15);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    // Minkowski metric(grid, m, a);
    SchwarzSchild metric(grid, m, a);
    // KerrSchild metric(grid, m, 0.9);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    size_t resX = 400;
    size_t resY = 300;
    size_t width = 26;
    size_t height = 19.5;
    Coord position(0,19,6);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 100 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Config:
    Config config =
    {
        .name = "ThinDiskE_" + metric.Name() + "_" + stencil.name + "_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z_"
              + std::to_string(sigma) + "_" + StreamingName(streamingType) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 50,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = false,
        .useCamera = true,
        .streamingType = streamingType,
        .initialDataType = InitialDataType::EandF,
    };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);
    radiation.sigma = sigma;

    // Initial Data:
    #pragma omp parallel for
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double radius = xyz.EuklNorm();
        double phi = xyz.Phi();

        // Disk:
        if(diskInner <= radius && radius <= diskOuter && abs(xyz[3]) < 0.9 * grid.dz)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.kappa0[ijk] = 0;
            radiation.kappa1[ijk] = 0;
            radiation.kappaA[ijk] = 0;
            radiation.eta[ijk] = 0;
            radiation.initialE[ijk] = 1;
        }
    }
    radiation.RunSimulation();
}

void ThinDisk2(size_t nx, size_t ny, size_t nz, Stencil stencil, int sigma, int simTime, StreamingType streamingType, string comment)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    Coord start(-15,-16,-10);
    Coord   end( 15, 22, 16);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    // Minkowski metric(grid, m, a);
    SchwarzSchild metric(grid, m, a);
    // KerrSchild metric(grid, m, 0.9);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    size_t resX = 400;
    size_t resY = 400;
    size_t width = 30;
    size_t height = 30;
    Coord position(0,19,1);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 95 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Config:
    Config config =
    {
        .name = "ThinDiskE2_" + metric.Name() + "_" + stencil.name + "_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z_"
              + std::to_string(sigma) + "_" + StreamingName(streamingType) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 50,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = false,
        .useCamera = true,
        .streamingType = streamingType,
        .initialDataType = InitialDataType::EandF,
    };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);
    radiation.sigma = sigma;

    // Initial Data:
    #pragma omp parallel for
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double radius = xyz.EuklNorm();
        double phi = xyz.Phi();

        // Disk:
        if(diskInner <= radius && radius <= diskOuter && abs(xyz[3]) < 0.9 * grid.dz)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.kappa0[ijk] = 0;
            radiation.kappa1[ijk] = 0;
            radiation.kappaA[ijk] = 0;
            radiation.eta[ijk] = 0;
            radiation.initialE[ijk] = 1;
        }
    }
    radiation.RunSimulation();
}
*/

int main(int argc, char *argv[])
{
    int n = 1;
    if (argc > 1)
        n = atoi(argv[1]);

    // SphereWaveAnalysis(n);       // Done
    // ShadowAnalysis(n);           // Done
    // BeamCrossingAnalysis(n);     // Done
    // DiffusionAnalysis(n);        // Done
    // MovingDiffusionAnalysis(n);  // Done
    // CurvedBeamAnalysis(n);       // Done
    // StencilAnalysis(n);          // Done

    // Thin Disk:
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,0),1,80, StreamingType::CurvedFixed, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,0),1,80, StreamingType::CurvedAdaptive, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,24),1,80, StreamingType::CurvedAdaptive, "");

    // High Res Runs:
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(31),1,80, StreamingType::CurvedFixed, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::CurvedFixed, "");
}