#include <iostream>
#include "../src/Radiation.h"
using namespace std;

void StraightBeamShadow(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 100 + 1 + 2;
    size_t ny = 50 + 1 + 2;
    size_t nz = 50 + 1 + 2;
    double dx = 2.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    double dz = 1.0 / (nz - 1.0 - 2.0);
    Coord start(0 - dx, -0.5 - dy, -0.5 - dz);
    Coord end(2 + dx, 0.5 + dy, 0.5 + dz);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    Camera camera;

    // Config:
    Config config =
        {
            .name = "Straight Beam Shadow 3d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 2,
            .writePeriod = 2, // write first and last frame only.
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
            .useCamera = false,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);

    // Initial Data:
    Coord center(1.0 / 3.0, 0.0, 0.0);
    double radius = 0.125;
    for (size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];
                radiation.initialKappa0[ijk] = 0;
                radiation.initialKappa1[ijk] = 0;
                radiation.initialKappaA[ijk] = 0;
                radiation.initialEta[ijk] = 0;
                if ((i <= 1) && (-0.25 < y && y < 0.25) && (-0.25 < z && z < 0.25))
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 1;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
                double dist = (center - xyz).EuklNorm();
                if (dist <= radius)
                    radiation.initialKappaA[ijk] = 1e10;
            }
    radiation.RunSimulation();
}

void Diffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl)
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 201;
    Coord start(-0.5, -0.5, -0.5);
    Coord end(0.5, 0.5, 0.5);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    Camera camera;

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (3.0 * kappa0) * (1.0 + 0.75 * PE); // Me
    // double D = 1.0 / (3.0 * kappa0) * (1.0 - 0.98 * PE) + grid.dx * grid.dx / grid.dt * 0.65; // Lukas
    // double D = 1.0 / (3.0 * kappa0); // analytical diffusion

    // Config:
    double t0 = 1;
    Config config =
        {
            // .name = "Diffusion 3d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string(nz) + "nz "  + Format(PE, 1, true) + "PE",
            .name = "Diffusion 3d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string(nz) + "nz " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE",
            .t0 = t0,
            .simTime = 1,
            .writePeriod = 1,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = true,
            .printToTerminal = true,
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
                radiation.initialKappa0[ijk] = kappa0;
                radiation.initialKappa1[ijk] = kappa1;
                radiation.initialKappaA[ijk] = 0;
                radiation.initialEta[ijk] = 0;
                radiation.isInitialGridPoint[ijk] = true;

                double E = pow(1.0 / t0, 3.0 / 2.0) * exp(-r * r / (4.0 * D * t0));
                radiation.initialE_LF[ijk] = E;
                radiation.initialFx_LF[ijk] = (x * E) / (2.0 * t0);
                radiation.initialFy_LF[ijk] = (y * E) / (2.0 * t0);
                radiation.initialFz_LF[ijk] = (z * E) / (2.0 * t0);
            }
    radiation.RunSimulation();
}

void CurvedBeam(Stencil stencil, StreamingType streamingType, double cfl, std::string comment)
{
    // Create Radiation object:
    size_t nx = 100 + 1 + 2;
    size_t ny = 80 + 1 + 2;
    size_t nz = 20 + 1 + 2;
    double dx = 5.0 / (nx - 1.0 - 2.0);
    double dy = 4.0 / (ny - 1.0 - 2.0);
    double dz = 1.0 / (nz - 1.0 - 2.0);
    Coord start(0 - dx, 0 - dy, -0.5 - dz);
    Coord end(5 + dx, 4 + dy, 0.5 + dz);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Camera:
    size_t resX = 100;
    size_t resY = 200;
    size_t width = 2;
    size_t height = 4;
    Coord position(2, 2, 0);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 0 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = -135 * degreeToRadians;
    glm::vec3 eulerAngles(angleX, angleY, angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Config:
    Config config =
        {
            .name = "Curved Beam 3d/" + metric.Name() + " " + stencil.name + " " + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" + Format(cfl, 2) + "cfl " + StreamingName(streamingType) + ((comment == "") ? "" : (" " + comment)),
            .simTime = 10.0,
            .writePeriod = 1,
            .updateSphericalHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
            .useCamera = false,
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
                radiation.initialKappa0[ijk] = 0;
                radiation.initialKappa1[ijk] = 0;
                radiation.initialKappaA[ijk] = 0;
                radiation.initialEta[ijk] = 0;
                if (x <= 1.1 * grid.dx && 3.00 <= y && y <= 3.50 && -0.25 <= z && z <= 0.25)
                {
                    Tensor4 uLF(1, 1, 0, 0);
                    uLF = NullNormalize(uLF, metric.GetMetric_ll(ijk));
                    Tensor3 vIF = Vec3ObservedByEulObs<LF, IF>(uLF, xyz, metric);

                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 1;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
            }
    radiation.RunSimulation();
}

/*
void SphereWave(Stencil stencil, StreamingType streamingType)
{
    // Grid, Metric, Stencil:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);

    // Config:
    Config config =
    {
        .name = "Sphere Wave 3d/" + StreamingName(streamingType) + "_" + stencil.name,
        .simTime = 0.7,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false,
        .streamingType = streamingType,
        .initialDataType = InitialDataType::EandF,
    };

    // Camera:
    size_t resX = 100;
    size_t resY = 100;
    size_t width = 2;
    size_t height = 2;
    Coord position(0,0.5,0.5);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 135 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);
    radiation.sigma = 1.0;

    // Initial Data:
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double r = xyz.EuklNorm();
        if (r < 0.1)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
        }
    }
    radiation.RunSimulation();
}


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
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
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
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
            radiation.initialE[ijk] = 1;
        }
    }
    radiation.RunSimulation();
}
*/

// Note:
// -TE changed from 1e-6 to 1e-4

// TODO:
// -fix kerr metric initial direction
// -try higher cfl condition
// -test Halo=1 vs Halo=2
// -UnstructuredMatrix use [] instead of .Row
int main(int argc, char *argv[])
{
    int n = 1;
    if (argc > 1)
        n = atoi(argv[1]);

    // StraightBeamShadow(LebedevStencil(35), StreamingType::FlatFixed, 0.75);
    // StraightBeamShadow(LebedevStencil(35, 4, 8, M_PI / 8.0), StreamingType::FlatAdaptive, 0.75);
    // StraightBeamShadow(LebedevStencil(23), StreamingType::FlatFixed, 0.75);
    // StraightBeamShadow(LebedevStencil(23, 4, 8, M_PI / 8.0), StreamingType::FlatAdaptive, 0.75);

    // Straight Beam Shadow:
    // if (n == 0)
    //    StraightBeamShadow(LebedevStencil(23), StreamingType::FlatFixed, 0.75);
    // if (n == 1)
    //    StraightBeamShadow(LebedevStencil(35), StreamingType::FlatFixed, 0.75);
    // if (n == 2)
    //    StraightBeamShadow(LebedevStencil(23), StreamingType::FlatAdaptive, 0.75);
    // if (n == 3)
    //    StraightBeamShadow(LebedevStencil(35), StreamingType::FlatAdaptive, 0.75);
    // if (n == 4)
    //    StraightBeamShadow(LebedevStencil(23, 4, 8, M_PI / 8.0), StreamingType::FlatAdaptive, 0.75);
    // if (n == 5)
    //    StraightBeamShadow(LebedevStencil(35, 4, 8, M_PI / 8.0), StreamingType::FlatAdaptive, 0.75);

    double lambda = 0;
    double cfl = 0.25;
    if (n == 0)
        Diffusion(LebedevStencil(21), StreamingType::FlatFixed, 100.0, lambda, cfl);
    // Diffusion(LebedevStencil(7, 2, 4, 0.3), StreamingType::FlatAdaptive, 100.0, lambda, cfl);
    if (n == 1)
        Diffusion(LebedevStencil(21), StreamingType::FlatFixed, 500.0, lambda, cfl);
    if (n == 2)
        Diffusion(LebedevStencil(21), StreamingType::FlatFixed, 1000.0, lambda, cfl);
    if (n == 3)
        Diffusion(LebedevStencil(21), StreamingType::FlatFixed, 10000.0, lambda, cfl);
    if (n == 4)
        Diffusion(LebedevStencil(21), StreamingType::FlatFixed, 100000.0, lambda, cfl);
    if (n == 5)
        Diffusion(LebedevStencil(21), StreamingType::FlatFixed, 1000000.0, lambda, cfl);

    // if (n == 0)
    //     CurvedBeam(LebedevStencil(23, 4, 8, M_PI / 8.0), StreamingType::CurvedFixed, 0.5, "500.1000 InterpGrid uint float");
    // if (n == 1)
    //     CurvedBeam(LebedevStencil(23, 4, 8, M_PI / 8.0), StreamingType::CurvedAdaptive, 0.5, "500.1000 InterpGrid uint float");

    // OLD:

    // Straight Beam:
    // StraightBeam(LebedevStencil(35),StreamingType::FlatFixed);
    // StraightBeam(LebedevStencil(35,5,6,M_PI/16.0),StreamingType::FlatAdaptive);
    // if (n==0) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.2);
    // if (n==1) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.3);
    // if (n==2) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.4);
    // if (n==3) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.5);
    // if (n==4) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.6);
    // if (n==5) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.7);
    // if (n==6) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.8);
    // if (n==7) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 0.9);
    // if (n==8) StraightBeam(LebedevStencil(23,5,8,M_PI/16.0),StreamingType::FlatAdaptive, 1.0);
    // if (n==1) StraightBeam(LebedevStencil(35,5,6,M_PI/16.0),StreamingType::FlatAdaptive);

    // Sphere Wave:
    // SphereWave(LebedevStencil(11), StreamingType::FlatFixed);
    // SphereWave(LebedevStencil(11), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(11,2,4,M_PI/8.0), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(11,3,8,M_PI/8.0), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(23), StreamingType::FlatFixed);
    // SphereWave(LebedevStencil(23), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(19,2,4,M_PI/8.0), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(19,3,8,M_PI/8.0), StreamingType::FlatAdaptive);

    // Curved Beam:
    // double dx = 0.050;
    // nx = deltaX / dx + 1
    // add 2 extra points for halo and extend deltaX by dx on both sides.
    // if (n==0) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(23)             , 20, 10, StreamingType::CurvedAdaptive, "");
    // if (n==1) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(23)             , 25, 10, StreamingType::CurvedAdaptive, "");
    // if (n==2) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(23)             , 30, 10, StreamingType::CurvedAdaptive, "");

    // if (n==0) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(23)              , 13, 10, StreamingType::CurvedFixed   , "");
    // if (n==1) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(23)              , 13, 10, StreamingType::CurvedAdaptive, "");
    // if (n==2) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(19,3,8,M_PI/ 8.0), 13, 10, StreamingType::CurvedAdaptive, "");
    // if (n==3) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(19,3,8,M_PI/16.0), 13, 10, StreamingType::CurvedAdaptive, "");
    // if (n==4) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(31)              , 22, 10, StreamingType::CurvedFixed   , "");
    // if (n==5) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(31)              , 22, 10, StreamingType::CurvedAdaptive, "");
    // if (n==6) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(29,3,8,M_PI/ 8.0), 22, 10, StreamingType::CurvedAdaptive, "");
    // if (n==7) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(29,3,8,M_PI/16.0), 22, 10, StreamingType::CurvedAdaptive, "");

    // if (n==0) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,4, 4,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==1) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,4,16,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==2) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,4,32,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==3) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,5, 4,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==4) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,5,16,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==5) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,5,32,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==6) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,6, 4,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==7) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,6,16,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");
    // if (n==8) CurvedBeam(dx, 5.0/dx+1, 4.0/dx+1, 1.0/dx+1, LebedevStencil(35,6,32,M_PI/16.0),30, 10, StreamingType::CurvedAdaptive   , "");

    // if (n==0) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35)              ,35, 10, StreamingType::CurvedFixed   , "");
    // if (n==1) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35)              ,35, 10, StreamingType::CurvedAdaptive, "");
    // if (n==2) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35,3,8,M_PI/ 8.0),35, 10, StreamingType::CurvedAdaptive, "");
    // if (n==3) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35,3,8,M_PI/16.0),35, 10, StreamingType::CurvedAdaptive, "");
    // if (n==4) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35)              ,40, 10, StreamingType::CurvedFixed   , "");
    // if (n==5) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35)              ,40, 10, StreamingType::CurvedAdaptive, "");
    // if (n==6) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35,3,8,M_PI/ 8.0),40, 10, StreamingType::CurvedAdaptive, "");
    // if (n==7) CurvedBeam(dx, 5.0/dx+1+2, 4.0/dx+1+2, 1.0/dx+1+2, LebedevStencil(35,3,8,M_PI/16.0),40, 10, StreamingType::CurvedAdaptive, "");

    // Thin Disk:
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,0),1,80, StreamingType::CurvedFixed, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,0),1,80, StreamingType::CurvedAdaptive, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,24),1,80, StreamingType::CurvedAdaptive, "");

    // High Res Runs:
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(31),1,80, StreamingType::CurvedFixed, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::CurvedFixed, "");

    // SphereWave(GaussLegendreStencil(5), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(15), StreamingType::FlatAdaptive);
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  35, 200, 10, StreamingType::CurvedAdaptive, "lebedevStreaming");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::CurvedAdaptive, "");
}