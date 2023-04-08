#include <iostream>
#include "../src/Radiation.h"
#include "../src/Radiation2.h"
#include "../src/Radiation3.h"
using namespace std;



void Test_Camera()
{
    // Black Hole with Thin Disk Camera Settings:
    {
        size_t resX = 400;
        size_t resY = 200;
        size_t width = 28;
        size_t height = 14;

        Coord position(0,20,3);
        double degreeToRadians = 2.0 * M_PI / 360.0;
        double angleX = 100 * degreeToRadians;
        double angleY = 0 * degreeToRadians;
        double angleZ = 0 * degreeToRadians;
        glm::vec3 eulerAngles(angleX,angleY,angleZ);

        Camera camera(resX, resY, width, height, position, eulerAngles);

        for(size_t ij=0; ij<camera.pixelCount; ij++)
        {
            size_t i = ij % camera.resX;
            size_t j = ij / camera.resX;
            camera.image[ij] = i / (resX-1.0) * j / (resY-1.0);
        }

        camera.WriteImagetoCsv(0, 0, "output");
    }
    
    // Curved Beam Close Settings:
    {
        size_t resX = 100;
        size_t resY = 50;
        size_t width = 2;
        size_t height = 1;

        Coord position(0,3/sqrt(2.0),3/sqrt(2.0));
        double degreeToRadians = 2.0 * M_PI / 360.0;
        double angleX = -145 * degreeToRadians;
        double angleY = 0 * degreeToRadians;
        double angleZ = 0 * degreeToRadians;
        glm::vec3 eulerAngles(angleX,angleY,angleZ);

        Camera camera(resX, resY, width, height, position, eulerAngles);

        for(size_t ij=0; ij<camera.pixelCount; ij++)
        {
            size_t i = ij % camera.resX;
            size_t j = ij / camera.resX;
            camera.image[ij] = i / (resX-1.0) * j / (resY-1.0);
        }

        camera.WriteImagetoCsv(0, 0, "output");
    }
}



void Test_Emission(size_t nx, size_t ny, size_t nz, size_t nTh, int sigma, int simTime)
{
    // Grid, Metric, Stencil:
    Coord start(-1,0,-0.1);
    Coord end(1,3.6,3.9);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);
    MyStencil stencil(nTh);
    LebedevStencil lebedevStencil(5);

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
    // Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedStatic);
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedDynamic);
    radiation.sigma = sigma;

    // Initial Data:
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
            radiation.isInitialGridPoint[ijk] = false;
            radiation.initialE[ijk] = 0;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 1;
        }
    }

    // Start simulation:
    Config config =
    {
        // .name = "Test_StreamCurvedStaticBeam",
        .name = "Emission " + metric.Name() + " " + std::to_string(stencil.nTh) + "." + std::to_string(stencil.nPh)
              + " s" + std::to_string(sigma) + " Leb" + std::to_string(lebedevStencil.nOrder) + " t" + std::to_string(simTime),
        .simTime = (double)simTime,
        .writeFrequency = 20,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}



void BlackHoleCsv()
{
    // Grid:
    Coord start(-12,-12,-12);
    Coord end(12,12,12);
    Grid grid(100,100,100,start,end);

    // Black Hole Geometry:
    double m = 1;
    double r = 2 * m;
    double diskInner = 3 * r;
    double diskOuter = 6 * r;

    ofstream file0("output/BlackHole.csv");
    ofstream file1("output/ThinDisk.csv");
    file0 << "#x,y,z,value\n";
    file1 << "#x,y,z,value\n";

    // Add one point with different value so min != max.
    file0 << "0,0,0,1" << "\n";
    file1 << "0,0,0,0" << "\n";

    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        Coord xyz = grid.xyz(i,j,k);
        double radius = xyz.EuklNorm();

        // Black Hole:
        if(radius <= r)
            file0 << xyz[1] << "," << xyz[2] << "," << xyz[3] << "," << 0 << "\n";
        // Disk:
        if(diskInner <= radius && radius <= diskOuter && abs(xyz[3]) < 0.2)
            file1 << xyz[1] << "," << xyz[2] << "," << xyz[3] << "," << 1 << "\n";
    }

    file0.close();
    file1.close();
}



void StreamFlatStaticSphereWave()
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    MyStencil stencil(7);
    LebedevStencil lebedevStencil(3);
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);
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
        .name = "Test_StreamFlatStaticSphereWave",
        .simTime = 1,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false
    };
    radiation.RunSimulation(config);
}



void StreamFlatDynamicSphereWave()
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    MyStencil stencil(7);
    LebedevStencil lebedevStencil(3);
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);
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
        .name = "Test_StreamFlatDynamicSphereWave",
        .simTime = 1,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false
    };
    radiation.RunSimulation(config);
}



void StreamFlatBeam()
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-2,0,0);
    Coord end(2,4,4);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    MyStencil stencil(15);
    LebedevStencil lebedevStencil(3);
    Camera camera;
    // Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);
    radiation.sigma = 100;

    // Initial Data:
    double beamWidth = 0.5;
    Coord beamCenter(0,3,0);
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double x = xyz[1];
        double y = xyz[2];
        double z = xyz[3];
        if (-0.5 < x && x < 0.5
          && 2.5 < y && y < 3.5
          && 0.2 < z && z < 0.3)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 1;
        }
    }

    // Start simulation:
    Config config =
    {
        // .name = "Test_StreamFlatDynamicBeam",
        .name = "Test_StreamFlatStaticBeam",
        .simTime = 4,
        .writeFrequency = 10,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false
    };
    radiation.RunSimulation(config);
}



void StreamFlatBeamCrossing()
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-2,-2,-2);
    Coord end(2,2,2);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    MyStencil stencil(15);
    LebedevStencil lebedevStencil(3);
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);
    // Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);
    radiation.sigma = 100;

    // Initial Data:
    double beamWidth = 0.5;
    Coord beamCenter(0,3,0);
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double x = xyz[1];
        double y = xyz[2];
        double z = xyz[3];
        if ( -0.5 < x && x <  0.5
          && -0.5 < y && y <  0.5
          && -2.0 < z && z < -1.8)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 1;
        }
        if (-0.5 < x && x <  0.5
          &&-2.0 < y && y < -1.8
          &&-0.5 < z && z <  0.5)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 1;
            radiation.initialNz[ijk] = 0;
        }
        if (-2.0 < x && x < -1.8
          &&-0.5 < y && y <  0.5
          &&-0.5 < z && z <  0.5)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 1;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
        }
    }

    // Start simulation:
    Config config =
    {
        // .name = "StreamFlatStaticBeamCrossing",
        .name = "StreamFlatDynamicBeamCrossing",
        .simTime = 4,
        .writeFrequency = 10,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false
    };
    radiation.RunSimulation(config);
}
void StreamCurvedBeamCrossing()
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(200-2,200-2,200-2);
    Coord end(200+2,200+2,200+2);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);
    MyStencil stencil(15);
    LebedevStencil lebedevStencil(5);

    // Camera:
    size_t resX = 100;
    size_t resY = 100;
    size_t width = 4;
    size_t height = 4;
    Coord position(200,200,201.5);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 180 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedDynamic);
    radiation.sigma = 100;

    // Initial Data:
    double beamWidth = 0.5;
    Coord beamCenter(0,3,0);
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double x = xyz[1];
        double y = xyz[2];
        double z = xyz[3];
        if ( 200-0.5 < x && x < 200+0.5
          && 200-0.5 < y && y < 200+0.5
          && 200-2.0 < z && z < 200-1.8)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 1;
        }
        if ( 200-0.5 < x && x < 200+0.5
          && 200-2.0 < y && y < 200-1.8
          && 200-0.5 < z && z < 200+0.5)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 1;
            radiation.initialNz[ijk] = 0;
        }
        if ( 200-2.0 < x && x < 200-1.8
          && 200-0.5 < y && y < 200+0.5
          && 200-0.5 < z && z < 200+0.5)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = 1;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
        }
    }

    // Start simulation:
    Config config =
    {
        // .name = "StreamFlatStaticBeamCrossing",
        .name = "StreamCurvedDynamicBeamCrossing",
        .simTime = 6,
        .writeFrequency = 20,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}



void StreamCurvedBeam(size_t nx, size_t ny, size_t nz, size_t nTh, int sigma, int simTime, StreamingType streamingType, std::string comment)
{
    // Grid, Metric, Stencil:
    Coord start(-1,0,-0.1);
    Coord end(1,3.6,3.9);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);   // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    MyStencil stencil(nTh);
    LebedevStencil lebedevStencil(5);

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
        .name = "Curved_Beam_" + metric.Name() + "_" + std::to_string(stencil.nTh) + "th" + std::to_string(stencil.nPh) + "ph_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" + std::to_string(simTime) + "t_"
              + "Leb" + std::to_string(lebedevStencil.nOrder) + "_" + StreamingName(streamingType) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 20,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}



enum IntensityProfile { Uniform, Linear, Squared, Cubic, SqFunc, CubFunc1, CubFunc2, CubFunc3, CubFunc4, CubFunc5, CubFunc6, CubFunc7, UniformToSquared1, UniformToSquared2, UniformToSquared3 };
std::string IntensityProfileName(int n)
{
	std::string name("unknown");
	switch (n)
	{
   		case  0: { name = "Uniform";    } break;
   		case  1: { name = "Linear";     } break;
   		case  2: { name = "Squared";    } break;
   		case  3: { name = "Cubic";      } break;
   		case  4: { name = "SqFunc";     } break;
   		case  5: { name = "CubFunc1";   } break;
   		case  6: { name = "CubFunc2";   } break;
   		case  7: { name = "CubFunc3";   } break;
   		case  8: { name = "CubFunc4";   } break;
   		case  9: { name = "CubFunc5";   } break;
   		case 10: { name = "CubFunc6";   } break;
   		case 11: { name = "CubFunc7";   } break;
   		case 12: { name = "UniformToSquared1"; } break;
   		case 13: { name = "UniformToSquared2"; } break;
   		case 14: { name = "UniformToSquared3"; } break;
		default: { ExitOnError("Invalid IntensityProfile"); }
	}
	return name;
}
void ThinDiskE(size_t nx, size_t ny, size_t nz, size_t nTh, int sigma, int simTime, StreamingType streamingType, IntensityProfile intensityProfile, string comment)
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
    Minkowski metric(grid, m, a);
    // SchwarzSchild metric(grid, m, a);
    MyStencil stencil(nTh);
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

            double alpha = metric.GetAlpha(xyz);
            double d = (diskOuter - xyz[2]) / (2.0 * diskOuter);
            if (intensityProfile == IntensityProfile::Uniform)
                radiation.initialE[ijk] = 1;
            if (intensityProfile == IntensityProfile::Linear)
                radiation.initialE[ijk] = 0.1 + 0.9*d;
            else if (intensityProfile == IntensityProfile::Squared)
                radiation.initialE[ijk] = 0.1 + 0.9*d*d;
            else if (intensityProfile == IntensityProfile::Cubic)
                radiation.initialE[ijk] = 0.1 + 0.9*d*d*d;

            else if (intensityProfile == IntensityProfile::SqFunc)
                radiation.initialE[ijk] = 1.7*d*d - 17.0/15.0*d + 0.2;
            else if (intensityProfile == IntensityProfile::CubFunc1)
                radiation.initialE[ijk] = 16*d*d*d - 12*d*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc2)
                radiation.initialE[ijk] = 24*d*d*d - 20*d*d + 2*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc3)
                radiation.initialE[ijk] = 32*d*d*d - 28*d*d + 4*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc4)
                radiation.initialE[ijk] = 28*d*d*d - 22*d*d+ 1*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc5)
                radiation.initialE[ijk] = 20*d*d*d - 6*d*d - 9*d + 4;
            else if (intensityProfile == IntensityProfile::CubFunc6)
                radiation.initialE[ijk] = 24*d*d*d - 14*d*d - 4*d + 3;
            else if (intensityProfile == IntensityProfile::CubFunc7)
                radiation.initialE[ijk] = 28*d*d*d - 20*d*d - 1*d + 2.5;

            else if (intensityProfile == IntensityProfile::UniformToSquared1)
            {
                if (d <= 0.5)
                    radiation.initialE[ijk] = 1;
                else
                    radiation.initialE[ijk] = 4.0*(d-0.5)*(d-0.5) + 1;
            }
            else if (intensityProfile == IntensityProfile::UniformToSquared2)
            {
                if (d <= 0.5)
                    radiation.initialE[ijk] = 1;
                else
                    radiation.initialE[ijk] = 8.0*(d-0.5)*(d-0.5) + 1;
            }
            else if (intensityProfile == IntensityProfile::UniformToSquared3)
            {
                if (d <= 0.5)
                    radiation.initialE[ijk] = 1;
                else
                    radiation.initialE[ijk] = 12.0*(d-0.5)*(d-0.5) + 1;
            }
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
        .name = "ThinDisk_E_" + metric.Name() + "_" + std::to_string(stencil.nTh) + "th" + std::to_string(stencil.nPh) + "ph_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" + std::to_string(simTime) + "t_"
              + "Leb" + std::to_string(lebedevStencil.nOrder) + "_" + StreamingName(streamingType) + "_" + IntensityProfileName(intensityProfile) + "_" + comment,
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
void ThinDiskEta(size_t nx, size_t ny, size_t nz, size_t nTh, int sigma, int simTime, StreamingType streamingType, IntensityProfile intensityProfile, string comment)
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
    Minkowski metric(grid, m, a);
    // SchwarzSchild metric(grid, m, a);
    MyStencil stencil(nTh);
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
            radiation.initialE[ijk] = 0;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;

            double alpha = metric.GetAlpha(xyz);
            double d = (diskOuter - xyz[2]) / (2.0 * diskOuter);
            if (intensityProfile == IntensityProfile::Uniform)
                radiation.initialEta[ijk] = 1 / alpha;
            if (intensityProfile == IntensityProfile::Linear)
                radiation.initialEta[ijk] = 0.1 + 0.9*d;
            else if (intensityProfile == IntensityProfile::Squared)
                radiation.initialEta[ijk] = 0.1 + 0.9*d*d;
            else if (intensityProfile == IntensityProfile::Cubic)
                radiation.initialEta[ijk] = 0.1 + 0.9*d*d*d;

            else if (intensityProfile == IntensityProfile::SqFunc)
                radiation.initialEta[ijk] = 1.7*d*d - 17.0/15.0*d + 0.2;
            else if (intensityProfile == IntensityProfile::CubFunc1)
                radiation.initialEta[ijk] = 16*d*d*d - 12*d*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc2)
                radiation.initialEta[ijk] = 24*d*d*d - 20*d*d + 2*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc3)
                radiation.initialEta[ijk] = 32*d*d*d - 28*d*d + 4*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc4)
                radiation.initialEta[ijk] = 28*d*d*d - 22*d*d+ 1*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc5)
                radiation.initialEta[ijk] = 20*d*d*d - 6*d*d - 9*d + 4;
            else if (intensityProfile == IntensityProfile::CubFunc6)
                radiation.initialEta[ijk] = 24*d*d*d - 14*d*d - 4*d + 3;
            else if (intensityProfile == IntensityProfile::CubFunc7)
                radiation.initialEta[ijk] = 28*d*d*d - 20*d*d - 1*d + 2.5;

            else if (intensityProfile == IntensityProfile::UniformToSquared1)
            {
                if (d <= 0.5)
                    radiation.initialEta[ijk] = 1;
                else
                    radiation.initialEta[ijk] = 4.0*(d-0.5)*(d-0.5) + 1;
            }
            else if (intensityProfile == IntensityProfile::UniformToSquared2)
            {
                if (d <= 0.5)
                    radiation.initialEta[ijk] = 1;
                else
                    radiation.initialEta[ijk] = 8.0*(d-0.5)*(d-0.5) + 1;
            }
            else if (intensityProfile == IntensityProfile::UniformToSquared3)
            {
                if (d <= 0.5)
                    radiation.initialEta[ijk] = 1;
                else
                    radiation.initialEta[ijk] = 12.0*(d-0.5)*(d-0.5) + 1;
            }
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
        .name = "ThinDisk_Eta_" + metric.Name() + "_" + std::to_string(stencil.nTh) + "th" + std::to_string(stencil.nPh) + "ph_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" + std::to_string(simTime) + "t_"
              + "Leb" + std::to_string(lebedevStencil.nOrder) + "_" + StreamingName(streamingType) + "_" + IntensityProfileName(intensityProfile) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 50,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = false,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}



void Test_Radiation2(int nOrder)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    Coord start(-1,-1,-1);
    Coord   end( 1, 1, 1);
    Grid grid(50, 50, 50, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, m, a);
    LebedevStencil intensityStencil(nOrder);
    LebedevStencil streamingStencil(5);

    // Camera:
    size_t resX = 200;
    size_t resY = 200;
    size_t width = 2;
    size_t height = 2;
    Coord position(0,0,0.9);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 180 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Radiation:
    Radiation2 radiation(metric, intensityStencil, streamingStencil, camera);
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
            radiation.initialNx[ijk] = xyz[1] / r;
            radiation.initialNy[ijk] = xyz[2] / r;
            radiation.initialNz[ijk] = xyz[3] / r;
        }
    }

    Config config =
    {
        .name = "Test_Radiation2(" + std::to_string(nOrder) + ")",
        .simTime = 1,
        .writeFrequency = 50,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = false,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}



void Test_Radiation3SphereWave(int nOrder, StreamingType streamingType)
{
    // Create Radiation object:
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    LebedevStencil stencil(nOrder);
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

    Radiation3 radiation(metric, stencil, lebedevStencil, camera, streamingType);
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
        .name = "Test_Radiation3SphereWave_" + StreamingName(streamingType) + "_" + to_string(nOrder),
        .simTime = 0.7,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = false
    };
    radiation.RunSimulation(config);
}
void Test_Radiation3CurvedBeam(size_t nx, size_t ny, size_t nz, int nOrder, int sigma, int simTime, StreamingType streamingType, std::string comment)
{
    // Grid, Metric, Stencil:
    Coord start(-0.6,1.0,-0.1);
    Coord end  ( 0.6,3.8, 3.5);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);   // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    LebedevStencil stencil(nOrder);
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
    Radiation3 radiation(metric, stencil, lebedevStencil, camera, streamingType);
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
        .name = "Test_Radiation3CurvedBeam_" + metric.Name() + "_" + std::to_string(stencil.nDir) +"nDir_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" + std::to_string(simTime) + "t_s"
              + std::to_string(sigma) + "_Leb" + std::to_string(lebedevStencil.nOrder) + "_" + StreamingName(streamingType) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 20,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = true
    };
    radiation.RunSimulation(config);
}

void Test_Radiation3ThinDiskE(size_t nx, size_t ny, size_t nz, size_t nOrder, int sigma, int simTime, StreamingType streamingType, IntensityProfile intensityProfile, string comment)
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
    LebedevStencil stencil(nOrder);
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
    Radiation3 radiation(metric, stencil, lebedevStencil, camera, streamingType);
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

            double alpha = metric.GetAlpha(xyz);
            double d = (diskOuter - xyz[2]) / (2.0 * diskOuter);
            if (intensityProfile == IntensityProfile::Uniform)
                radiation.initialE[ijk] = 1;
            if (intensityProfile == IntensityProfile::Linear)
                radiation.initialE[ijk] = 0.1 + 0.9*d;
            else if (intensityProfile == IntensityProfile::Squared)
                radiation.initialE[ijk] = 0.1 + 0.9*d*d;
            else if (intensityProfile == IntensityProfile::Cubic)
                radiation.initialE[ijk] = 0.1 + 0.9*d*d*d;

            else if (intensityProfile == IntensityProfile::SqFunc)
                radiation.initialE[ijk] = 1.7*d*d - 17.0/15.0*d + 0.2;
            else if (intensityProfile == IntensityProfile::CubFunc1)
                radiation.initialE[ijk] = 16*d*d*d - 12*d*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc2)
                radiation.initialE[ijk] = 24*d*d*d - 20*d*d + 2*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc3)
                radiation.initialE[ijk] = 32*d*d*d - 28*d*d + 4*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc4)
                radiation.initialE[ijk] = 28*d*d*d - 22*d*d+ 1*d + 2;
            else if (intensityProfile == IntensityProfile::CubFunc5)
                radiation.initialE[ijk] = 20*d*d*d - 6*d*d - 9*d + 4;
            else if (intensityProfile == IntensityProfile::CubFunc6)
                radiation.initialE[ijk] = 24*d*d*d - 14*d*d - 4*d + 3;
            else if (intensityProfile == IntensityProfile::CubFunc7)
                radiation.initialE[ijk] = 28*d*d*d - 20*d*d - 1*d + 2.5;

            else if (intensityProfile == IntensityProfile::UniformToSquared1)
            {
                if (d <= 0.5)
                    radiation.initialE[ijk] = 1;
                else
                    radiation.initialE[ijk] = 4.0*(d-0.5)*(d-0.5) + 1;
            }
            else if (intensityProfile == IntensityProfile::UniformToSquared2)
            {
                if (d <= 0.5)
                    radiation.initialE[ijk] = 1;
                else
                    radiation.initialE[ijk] = 8.0*(d-0.5)*(d-0.5) + 1;
            }
            else if (intensityProfile == IntensityProfile::UniformToSquared3)
            {
                if (d <= 0.5)
                    radiation.initialE[ijk] = 1;
                else
                    radiation.initialE[ijk] = 12.0*(d-0.5)*(d-0.5) + 1;
            }
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
        .name = std::to_string(stencil.nDir) +"nDir_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" + std::to_string(simTime) + "t_s"
              + std::to_string(sigma) + "_Leb" + std::to_string(lebedevStencil.nOrder) + "_" + StreamingName(streamingType) + "_"
              + IntensityProfileName(intensityProfile) + "_" + comment,
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
// -fix bh rendering with dynamic stencil. Camera broken?
// -test ComputeMomentsIF with velocty rotation vs moment rotation.
int main(int argc, char *argv[])
{
    int n;
    if(argc > 1)
        n = atoi(argv[1]);

    // Tests:
    // Test_Camera();
    // Test_Emission( 25, 45, 50,  20,    15,10);
    // Test_Radiation2(7);
    // Test_Radiation2(9);
    // Test_Radiation2(11);
    // Test_Radiation2(13);
    // Test_Radiation3SphereWave(9, StreamingType::FlatStatic);
    Test_Radiation3SphereWave(9, StreamingType::FlatDynamic);
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  35, 200, 10, StreamingType::CurvedDynamic, "lebedevStreaming");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "");
    
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  17, 1,80, StreamingType::CurvedStatic  ,IntensityProfile::Uniform, "Ground0");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  17, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "BruteForce");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  17, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "minMaxDot");

    
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  35, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "NearestDirection");

    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  29, 1,80, StreamingType::CurvedStatic  ,IntensityProfile::Uniform, "BenchMarkStatic");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  7, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "BenchMarkDynamic");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  35, 1,80, StreamingType::CurvedStatic  ,IntensityProfile::Uniform, "BenchMarkStatic");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  35, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "BenchMarkDynamic");


    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  23, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "minMaxDot");
    
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  35, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "");
    // Visualization data:
    // BlackHoleCsv();

    // Runs:
    // StreamFlatBeam();
    // StreamFlatBeamCrossing();
    // StreamCurvedBeamCrossing();
    // StreamFlatStaticSphereWave();
    // StreamFlatDynamicSphereWave();

  //StreamCurvedBeam( nx, ny, nz, nTh, sigma, simTime, streamingType, comment)

    // StreamCurvedBeam( 26, 46, 51,  15,    15,10, StreamingType::CurvedStatic, "InvDistInterp_limited");
    // StreamCurvedBeam( 26, 46, 51,  15,    15,10, StreamingType::CurvedDynamic, "InvDistInterp_limited");
    // StreamCurvedBeam( 26, 46, 51,  19,   20,10);
    // StreamCurvedBeam( 51, 91, 101,  15,    15,10);
    // StreamCurvedBeam( 51, 91, 101,  19,    20,10);
    

 // ThinDiskE( nx, ny, nz, nTh, sigma, simTime, streamingType, intensityProfile, comment)
    // n€[3,10] 8 takes about 13.5h
    // ThinDiskE(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::FlatStatic  ,IntensityProfile::Uniform, "Bicubic");
    // ThinDiskE(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::FlatDynamic ,IntensityProfile::Uniform, "Bicubic");
    // ThinDiskEta(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::FlatDynamic ,IntensityProfile::Uniform, "Bicubic");
    // ThinDiskE(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::CurvedDynamic ,IntensityProfile::Uniform,           "harmonicRefactor");
    // ThinDiskE(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::CurvedStatic  ,IntensityProfile::Uniform,           "harmonicRefactor");
}