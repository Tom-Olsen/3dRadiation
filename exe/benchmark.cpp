#include <iostream>
#include "../src/Radiation.h"
using namespace std;

// Macros:
#define WRITE_DATA false
#define PRINT_SETUP false
#define PRINT_PROGRESS false
#define PRINT_RESULTS false
#define SAVE_ID false

double FlatTest(LebedevStencil stencil, StreamingType streamingType, double cfl, double kappa0, double kappa1, double kappaA, double eta, double ux)
{
    // Create Radiation object:
    size_t nx = 51;
    size_t ny = 51;
    size_t nz = 51;
    Coord start(4, 4, 4);
    Coord end(8, 8, 8);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Config:
    Config config =
        {
            .name = "FlatTest/" + StreamingName(streamingType) + " " + stencil.name + (kappa0 > 0 ? " with collisions" : " no collisions"),
            .t0 = 0,
            .simTime = 50 * grid.dt,
            .writePeriod = 100 * grid.dt,
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
    Coord center = (start + end) / 2;
    for(size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double r = (xyz - center).EuklNorm();
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];

                radiation.kappa0[ijk] = kappa0;
                radiation.kappa1[ijk] = kappa1;
                radiation.ux[ijk] = ux;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;

                if (0.5 < r && r < 1.0)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.eta[ijk] = eta;
                    radiation.kappaA[ijk] = kappaA;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
            }
    radiation.RunSimulation();

    std::cout << "average lambda itt count: " << radiation.logger.averageItterationCount << std::endl;
    int timeSteps = radiation.logger.timeSteps;
    double time = radiation.logger.ComputationTime();
    int nxyz = radiation.logger.metric.grid.nxyz;
    return nxyz * timeSteps / (1e6 * time);
}
double CurvedTest(LebedevStencil stencil, StreamingType streamingType, double cfl, double kappa0, double kappa1, double kappaA, double eta, double ux)
{
    // Create Radiation object:
    size_t nx = 51;
    size_t ny = 51;
    size_t nz = 51;
    Coord start(4, 4, 4);
    Coord end(8, 8, 8);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0);
    LebedevStencil streamingStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;

    // Config:
    Config config =
        {
            .name = "CurvedTest/" + StreamingName(streamingType) + " " + stencil.name + (kappa0 > 0 ? " with collisions" : " no collisions"),
            .t0 = 0,
            .simTime = 50 * grid.dt,
            .writePeriod = 100 * grid.dt,
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
    Coord center = (start + end) / 2;
    for(size_t k = 0; k < grid.nz; k++)
        for (size_t j = 0; j < grid.ny; j++)
            for (size_t i = 0; i < grid.nx; i++)
            {
                size_t ijk = grid.Index(i, j, k);
                Coord xyz = grid.xyz(i, j, k);
                double r = (xyz - center).EuklNorm();
                double x = xyz[1];
                double y = xyz[2];
                double z = xyz[3];

                radiation.kappa0[ijk] = kappa0;
                radiation.kappa1[ijk] = kappa1;
                radiation.ux[ijk] = ux;
                radiation.uy[ijk] = 0;
                radiation.uz[ijk] = 0;

                if (0.5 < r && r < 1.0)
                {
                    radiation.isInitialGridPoint[ijk] = true;
                    radiation.eta[ijk] = eta;
                    radiation.kappaA[ijk] = kappaA;
                    radiation.initialE_LF[ijk] = 1;
                    radiation.initialFx_LF[ijk] = 0;
                    radiation.initialFy_LF[ijk] = 0;
                    radiation.initialFz_LF[ijk] = 0;
                }
            }
    radiation.RunSimulation();

    std::cout << "average lambda itt count: " << radiation.logger.averageItterationCount << std::endl;
    int timeSteps = radiation.logger.timeSteps;
    double time = radiation.logger.ComputationTime();
    int nxyz = radiation.logger.metric.grid.nxyz;
    return nxyz * timeSteps / (1e6 * time);
}



void Benchmark(StreamingType streamingType, bool collisions)
{
    // Test Config:
    bool fixed    = streamingType == StreamingType::FlatFixed    || streamingType == StreamingType::CurvedFixed;
    bool adaptive = streamingType == StreamingType::FlatAdaptive || streamingType == StreamingType::CurvedAdaptive;
    bool flat     = streamingType == StreamingType::FlatFixed    || streamingType == StreamingType::FlatAdaptive;
    bool curved   = streamingType == StreamingType::CurvedFixed  || streamingType == StreamingType::CurvedAdaptive;

    // Parameters:
    double cfl = 0.9;
    int N = 10;
    int numberOfOrders = 14;
    int order[numberOfOrders] = {9, 11, 15, 17, 19, 21, 23, 29, 31, 35, 41, 47, 53, 59};
    Tensor3 results[numberOfOrders];
    double lups;

    // Collision Parameters:
    double kappa0 = 100;
    double kappa1 = 20;
    double kappaA = 10;
    double eta = 10;
    double ux = 0.5;

    // Output name of the test:
    string name;
    if (  flat &&    fixed && !collisions) name = "Flat Fixed no collisions";
    if (  flat &&    fixed &&  collisions) name = "Flat Fixed with collisions";
    if (  flat && adaptive && !collisions) name = "Flat Adaptive no collisions";
    if (  flat && adaptive &&  collisions) name = "Flat Adaptive with collisions";
    if (curved &&    fixed && !collisions) name = "Curved Fixed no collisions";
    if (curved &&    fixed &&  collisions) name = "Curved Fixed with collisions";
    if (curved && adaptive && !collisions) name = "Curved Adaptive no collisions";
    if (curved && adaptive &&  collisions) name = "Curved Adaptive with collisions";
    std::cout << name << std::endl;
    for(int i=0; i<numberOfOrders; i++)
    {
        double maxLups = -1;
        LebedevStencil stencil = fixed ? LebedevStencil(order[i], 0.00, 0.00) : LebedevStencil(order[i], 0.15, 0.00);
        std::cout << "name = " << stencil.name << std::endl;
        std::cout << "order = " << stencil.nOrder << std::endl;
        std::cout << "ndir = " << stencil.nDir << std::endl;
        std::cout << "nGhost = " << stencil.nGhost << std::endl;
        for (int j = 0; j < N; j++)
        {
            if (  flat &&    fixed && !collisions) lups =   FlatTest(stencil, StreamingType::FlatFixed     , cfl, 0, 0, 0, 0, 0);
            if (  flat &&    fixed &&  collisions) lups =   FlatTest(stencil, StreamingType::FlatFixed     , cfl, kappa0, kappa1, kappaA, eta, ux);
            if (  flat && adaptive && !collisions) lups =   FlatTest(stencil, StreamingType::FlatAdaptive  , cfl, 0, 0, 0, 0, 0);
            if (  flat && adaptive &&  collisions) lups =   FlatTest(stencil, StreamingType::FlatAdaptive  , cfl, kappa0, kappa1, kappaA, eta, ux);
            if (curved &&    fixed && !collisions) lups = CurvedTest(stencil, StreamingType::CurvedFixed   , cfl, 0, 0, 0, 0, 0);
            if (curved &&    fixed &&  collisions) lups = CurvedTest(stencil, StreamingType::CurvedFixed   , cfl, kappa0, kappa1, kappaA, eta, ux);
            if (curved && adaptive && !collisions) lups = CurvedTest(stencil, StreamingType::CurvedAdaptive, cfl, 0, 0, 0, 0, 0);
            if (curved && adaptive &&  collisions) lups = CurvedTest(stencil, StreamingType::CurvedAdaptive, cfl, kappa0, kappa1, kappaA, eta, ux);

            std::cout << "Benchmark: " << lups << "MLUPS" << std::endl;
            maxLups = std::max(maxLups, lups);
        }
        results[i] = Tensor3(stencil.nDir, maxLups, 0);
        std::cout << "max: " << maxLups << std::endl;
        std::cout << std::endl;
    }
    std::cout << "#" << name << std::endl;
    std::cout << "#order, MLUPS: " << std::endl;
    for(int i=0; i<numberOfOrders; i++)
        std::cout << results[i][1] << ", " << results[i][2] << std::endl;
}



int main(int argc, char *argv[])
{
    int n = 1;
    if (argc > 1)
        n = atoi(argv[1]);

    double cfl = 0.9;
    double kappa0 = 100;
    double kappa1 = 20;
    double kappaA = 10;
    double eta = 10;
    double ux = 0.5;

    // if (n == 0) Benchmark(StreamingType::FlatFixed     , false);
    // if (n == 1) Benchmark(StreamingType::FlatFixed     , true);
    // if (n == 2) Benchmark(StreamingType::FlatAdaptive  , false);
    // if (n == 3) Benchmark(StreamingType::FlatAdaptive  , true);
    // if (n == 4) Benchmark(StreamingType::CurvedFixed   , false);
    // if (n == 5) Benchmark(StreamingType::CurvedFixed   , true);
    // if (n == 6) Benchmark(StreamingType::CurvedAdaptive, false);
    // if (n == 7) Benchmark(StreamingType::CurvedAdaptive, true);

    // Scaling analysis. Run this with different OMP_NUM_THREADS:
    // OMP_NUM_THREADS=N ./benchmark.out
    // const char* omp_num_threads_env = std::getenv("OMP_NUM_THREADS");
    // int omp_num_threads = 1;
    // if (omp_num_threads_env != nullptr)
    // {
        // omp_num_threads = std::stoi(omp_num_threads_env);
        // omp_set_num_threads(omp_num_threads);
    // }
    // double maxLups = -1;
    // for(int i=0; i<10; i++)
    // {
        // double lups =  CurvedTest(LebedevStencil(19, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl, kappa0, kappa1, kappaA, eta, ux);
        // maxLups = std::max(maxLups, lups);
    // }
    // std::cout << "# Threads, MLUPS" << std::endl;
    // std::cout << omp_num_threads_env << ", " << maxLups << std::endl;
    
    // if (n == 0)   FlatTest(LebedevStencil(19, 0.00, 0.00, 0.00), StreamingType::FlatFixed     , cfl, 0, 0, 0, 0, 0);
    // if (n == 1)   FlatTest(LebedevStencil(19, 0.00, 0.00, 0.00), StreamingType::FlatFixed     , cfl, kappa0, kappa1, kappaA, eta, ux);
    // if (n == 2)   FlatTest(LebedevStencil(19, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive  , cfl, 0, 0, 0, 0, 0);
    // if (n == 3)   FlatTest(LebedevStencil(19, 0.15, 0.00, 0.00), StreamingType::FlatAdaptive  , cfl, kappa0, kappa1, kappaA, eta, ux);
    // if (n == 4) CurvedTest(LebedevStencil(19, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, 0, 0, 0, 0, 0);
    // if (n == 5) CurvedTest(LebedevStencil(19, 0.00, 0.00, 0.00), StreamingType::CurvedFixed   , cfl, kappa0, kappa1, kappaA, eta, ux);
    // if (n == 6) CurvedTest(LebedevStencil(19, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl, 0, 0, 0, 0, 0);
    // if (n == 7) CurvedTest(LebedevStencil(19, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl, kappa0, kappa1, kappaA, eta, ux);

    // Single run (perf):
    CurvedTest(LebedevStencil(19, 0.15, 0.00, 0.00), StreamingType::CurvedAdaptive, cfl, kappa0, kappa1, kappaA, eta, ux);
}