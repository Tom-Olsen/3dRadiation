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

    ofstream file0("output/BlackHole.txt");
    ofstream file1("output/ThinDisk.txt");
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



void Test_Radiation3SphereWave(Stencil stencil, StreamingType streamingType)
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
        .name = "Test_Radiation3SphereWave_" + StreamingName(streamingType) + "_" + stencil.name,
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
void Test_Radiation3CurvedBeam(size_t nx, size_t ny, size_t nz, Stencil stencil, int sigma, int simTime, StreamingType streamingType, std::string comment)
{
    // Grid, Metric, Stencil:
    Coord start(-0.6,1.0,-0.1);
    Coord end  ( 0.6,3.8, 3.5);
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
        .name = "Radiation3CurvedBeam_" + metric.Name() + "_" + stencil.name + "_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z_"
              + std::to_string(sigma) + "_" + StreamingName(streamingType) + "_" + comment,
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

void Test_Radiation3ThinDiskE(size_t nx, size_t ny, size_t nz, Stencil stencil, int sigma, int simTime, StreamingType streamingType, IntensityProfile intensityProfile, string comment)
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
        .name = "Radiation3ThinDiskE_" + metric.Name() + "_" + stencil.name + "_"
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

void Test_Radiation3ThinDiskEotherCamera(size_t nx, size_t ny, size_t nz, Stencil stencil, int sigma, int simTime, StreamingType streamingType, string comment)
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
        .name = "Test_Radiation3ThinDiskE2_" + metric.Name() + "_" + stencil.name + "_"
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
// -fix bh rendering with dynamic stencil.
int main(int argc, char *argv[])
{
    int n = 1;
    if(argc > 1)
        n = atoi(argv[1]);
        
    // Sphere Wave:
    // Test_Radiation3SphereWave(GaussLegendreStencil( 5), StreamingType::FlatStatic);
    // Test_Radiation3SphereWave(GaussLegendreStencil(15), StreamingType::FlatStatic);
    // Test_Radiation3SphereWave(GaussLegendreStencil( 5), StreamingType::FlatDynamic);
    // Test_Radiation3SphereWave(GaussLegendreStencil(15), StreamingType::FlatDynamic);
    // Test_Radiation3SphereWave(LebedevStencil(11), StreamingType::FlatStatic);
    // Test_Radiation3SphereWave(LebedevStencil(35), StreamingType::FlatStatic);
    // Test_Radiation3SphereWave(LebedevStencil(11), StreamingType::FlatDynamic);
    // Test_Radiation3SphereWave(LebedevStencil(35), StreamingType::FlatDynamic);

    // Thin Disk:
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil( 5),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil( 7),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil( 9),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil(15),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil( 5),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil( 7),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil( 9),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,GaussLegendreStencil(15),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(11),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(17),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(21),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(11),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(17),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(21),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::CurvedDynamic, IntensityProfile::Uniform, "");

    // Tests:
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::FlatStatic, IntensityProfile::Uniform, "flat");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::FlatDynamic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskEotherCamera(n*15+1,n*19+1,n*13+1,LebedevStencil(35),1,80, StreamingType::CurvedStatic, "");

    // Curved Beam:
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(35),  50, 10, StreamingType::CurvedStatic, "");
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(35), 100, 10, StreamingType::CurvedStatic, "");
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(35), 200, 10, StreamingType::CurvedStatic, "");
    
    Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(35),  50, 10, StreamingType::CurvedDynamic, "");
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(35), 100, 10, StreamingType::CurvedDynamic, "");
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  LebedevStencil(35), 200, 10, StreamingType::CurvedDynamic, "");


    // High Res Runs:
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(31),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,LebedevStencil(35),1,80, StreamingType::CurvedStatic, IntensityProfile::Uniform, "");

    // Test_Radiation3SphereWave(GaussLegendreStencil(5), StreamingType::FlatDynamic);
    // Test_Radiation3SphereWave(LebedevStencil(15), StreamingType::FlatDynamic);
    // Test_Radiation3CurvedBeam(15*n+1,35*n+1,45*n+1,  35, 200, 10, StreamingType::CurvedDynamic, "lebedevStreaming");
    // Test_Radiation3ThinDiskE(n*21+1,n*27+1,n*12+1,  9, 1,80, StreamingType::CurvedDynamic  ,IntensityProfile::Uniform, "");
}