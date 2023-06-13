#include <iostream>
#include "../src/Radiation.h"
using namespace std;



void SphereWave(Stencil stencil, StreamingType streamingType)
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil lebedevStencil(5);
    
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

    Radiation radiation(metric, stencil, lebedevStencil, camera, streamingType);
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

    // Start simulation:
    Config config =
    {
        .name = "SphereWave_" + StreamingName(streamingType) + "_" + stencil.name,
        .simTime = 0.7,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}
void CurvedBeam(size_t nx, size_t ny, size_t nz, Stencil stencil, int sigma, int simTime, StreamingType streamingType, std::string comment)
{
    // Grid, Metric, Stencil:
    Coord start(-0.6, 1.0, -0.1);
    Coord end  ( 0.6, 3.8,  3.5);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);   // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    LebedevStencil lebedevStencil(5);

    // stencil.connectedTriangles.Print();
    // return;

    // Camera:
    size_t resX = 100;
    size_t resY = 200;
    size_t width = 2;
    size_t height = 4;
    Coord position(0,2,2);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = -135 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, camera, streamingType);
    radiation.sigma = sigma;

    // Initial Data:
    PARALLEL_FOR(3)
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double x = xyz[1];
        double y = xyz[2];
        double z = xyz[3];
        if (-0.25 < x && x < 0.25
          && 3.00 < y && y < 3.50
          && z < 0.00)
        {
            Tensor4 uLF(1,0,0,1);
            uLF = NullNormalize(uLF,metric.GetMetric_ll(ijk));
            Tensor3 vIF = Vec3ObservedByEulObs<LF,IF>(uLF, xyz, metric);

            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = vIF[1];
            radiation.initialNy[ijk] = vIF[2];
            radiation.initialNz[ijk] = vIF[3];
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
        }
    }

    // Start simulation:
    Config config =
    {
        .name = "CurvedBeam_" + metric.Name() + "_" + stencil.name + "_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z_"
              + std::to_string(sigma) + "s_" + StreamingName(streamingType) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 20,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false
    };
    radiation.RunSimulation(config);
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

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, camera, streamingType);
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
    

    // Get current time and date:
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%H.%M.%S - %d.%m.%Y");
	string date = oss.str();

    // Start simulation:
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
        .useCamera = true
    };
    radiation.RunSimulation(config);
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

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, camera, streamingType);
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
    

    // Get current time and date:
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%H.%M.%S - %d.%m.%Y");
	string date = oss.str();

    // Start simulation:
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
        .useCamera = true
    };
    radiation.RunSimulation(config);
}



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
    if(argc > 1)
        n = atoi(argv[1]);
        
    // Sphere Wave:
    // SphereWave(GaussLegendreStencil( 5), StreamingType::FlatFixed);
    // SphereWave(GaussLegendreStencil(15), StreamingType::FlatFixed);
    // SphereWave(GaussLegendreStencil( 5), StreamingType::FlatAdaptive);
    // SphereWave(GaussLegendreStencil(15), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(11), StreamingType::FlatFixed);
    // SphereWave(LebedevStencil(35), StreamingType::FlatFixed);
    // SphereWave(LebedevStencil(11), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(35), StreamingType::FlatAdaptive);

    // Thin Disk:
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,0),1,80, StreamingType::CurvedFixed, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,0),1,80, StreamingType::CurvedAdaptive, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(23,24),1,80, StreamingType::CurvedAdaptive, "");

    // Curved Beam:
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0), 30, 10, StreamingType::CurvedFixed, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0), 40, 10, StreamingType::CurvedFixed, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0), 50, 10, StreamingType::CurvedFixed, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0), 60, 10, StreamingType::CurvedFixed, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0), 70, 10, StreamingType::CurvedFixed, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0),  30, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0),  40, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0),  50, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0),  60, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,0),  70, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,24),  10, 10, StreamingType::CurvedAdaptive, "Halo2");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,24),  30, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,24),  40, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,24),  50, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,24),  60, 10, StreamingType::CurvedAdaptive, "");
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(23,24),  70, 10, StreamingType::CurvedAdaptive, "");


    // High Res Runs:
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(31),1,80, StreamingType::CurvedFixed, "");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::CurvedFixed, "");

    // SphereWave(GaussLegendreStencil(5), StreamingType::FlatAdaptive);
    // SphereWave(LebedevStencil(15), StreamingType::FlatAdaptive);
    // CurvedBeam(15*n+1,35*n+1,45*n+1,  35, 200, 10, StreamingType::CurvedAdaptive, "lebedevStreaming");
    // ThinDisk1(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::CurvedAdaptive, "");
}