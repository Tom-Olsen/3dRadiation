#include <fstream>
#include "../src/Stencil.h"
#include "../src/GeodesicEquationSolver.h"
#include "../src/Camera.h"



void BlackHoleCsv(Metric& metric)
{
    std::ofstream fileOut(OUTPUTDIR + "BlackHole.csv");

    int nx = 50;
    int ny = 50;
    int nz = 50;
    double r = 2.0 * metric.m;
    Coord start(-r,-r,-r);
    Coord end(r,r,r);
    Grid grid(nx,ny,nz,start,end);

    // Config:
    fileOut << "#nx=" << grid.nx << "\n";
    fileOut << "#ny=" << grid.ny << "\n";
    fileOut << "#nz=" << grid.nz << "\n";
    fileOut << "#startx=" << grid.startx << "\n";
    fileOut << "#starty=" << grid.starty << "\n";
    fileOut << "#startz=" << grid.startz << "\n";
    fileOut << "#endx=" << grid.endx << "\n";
    fileOut << "#endy=" << grid.endy << "\n";
    fileOut << "#endz=" << grid.endz << "\n";

    // Data:
    fileOut << "#x,y,z,value\n";
	for(size_t k = 0; k < grid.nz; k++)
	for(size_t j = 0; j < grid.ny; j++)
	for(size_t i = 0; i < grid.nx; i++)
    {
        Coord xyz = grid.xyz(i,j,k);
        double radius = xyz.EuklNorm();
        fileOut << xyz[1] << "," << xyz[2] << "," << xyz[3] << ",";
        if(radius <= r)
            fileOut << 1.0 << "\n";
        else
            fileOut << 0.0 << "\n";
    }
    fileOut.close();
}



void ThinDiskCsv(Metric& metric, double diskInner, double diskOuter)
{
    std::ofstream fileOut(OUTPUTDIR + "ThinDisk.csv");

    int nx = 100;
    int ny = 100;
    int nz = 10;
    double r = 2.0 * metric.m;
    Coord start(-diskOuter,-diskOuter,-0.2);
    Coord end(diskOuter,diskOuter,0.2);
    Grid grid(nx,ny,nz,start,end);

    // Config:
    fileOut << "#nx=" << grid.nx << "\n";
    fileOut << "#ny=" << grid.ny << "\n";
    fileOut << "#nz=" << grid.nz << "\n";
    fileOut << "#startx=" << grid.startx << "\n";
    fileOut << "#starty=" << grid.starty << "\n";
    fileOut << "#startz=" << grid.startz << "\n";
    fileOut << "#endx=" << grid.endx << "\n";
    fileOut << "#endy=" << grid.endy << "\n";
    fileOut << "#endz=" << grid.endz << "\n";

    // Data:
    fileOut << "#x,y,z,value\n";
	for(size_t k = 0; k < grid.nz; k++)
	for(size_t j = 0; j < grid.ny; j++)
	for(size_t i = 0; i < grid.nx; i++)
    {
        Coord xyz = grid.xyz(i,j,k);
        double radius = xyz.EuklNorm();
        fileOut << xyz[1] << "," << xyz[2] << "," << xyz[3] << ",";
        if(diskInner <= radius && radius <= diskOuter && abs(xyz[3]) < 0.1)
            fileOut << 1.0 << "\n";
        else
            fileOut << 0.0 << "\n";
    }
    fileOut.close();
}



void RaysFromCamera(Metric& metric, Camera& camera, double diskInner, double diskOuter, int steps)
{
    RealBuffer image; image.resize(camera.pixelCount);
    RealBuffer pixelX; pixelX.resize(camera.pixelCount);
    RealBuffer pixelY; pixelY.resize(camera.pixelCount);
    RealBuffer pixelZ; pixelZ.resize(camera.pixelCount);
    double delta = 0.1 * diskInner;

    PARALLEL_FOR(2)
    for(int j=0; j<camera.resY; j++)
    for(int i=0; i<camera.resX; i++)
    {
        // Emission point:
        Coord xyz0 = camera.xyz(i,j);
        int index = camera.Index(i,j);

        // Initial data for photon:
        Coord xyz = xyz0;
		double alpha = metric.GetAlpha(xyz);
        double s = 1;
        Tensor3 uLF3 = -camera.lookDirection;
        Tensor4 uLF(1,uLF3[1],uLF3[2],uLF3[3]);
        uLF = NullNormalize(uLF,metric.GetMetric_ll(xyz));
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF, xyz, metric);

		// Solve geodesic equation backwards:
        bool trackDiskHit = true;
        for(int k=0; k<steps; k++)
        {
            // Break conditions:
		    if(metric.InsideBH(xyz))
                break;
            if(metric.grid.OutsideDomain(xyz))
                break;
            
            // Check if disk edge is hit:
            double radius = xyz.EuklNorm();
            if(abs(xyz[3]) < 0.1 && trackDiskHit)
            {
                if(diskInner < radius && radius < diskInner + delta)
                {
                    image[index] = 1.0;
                    pixelX[index] = xyz0[1];
                    pixelY[index] = xyz0[2];
                    pixelZ[index] = xyz0[3];
                    trackDiskHit = false;
                    break;
                }
                if(diskOuter - delta < radius && radius < diskOuter)
                {
                    image[index] = 1.0;
                    pixelX[index] = xyz0[1];
                    pixelY[index] = xyz0[2];
                    pixelZ[index] = xyz0[3];
                    trackDiskHit = false;
                    break;
                }
            }

    	    s *= RK45_GeodesicEquation<-1>(0.5 * metric.grid.dt, xyz, vLF, metric);
        }
        // Pixel is 1.0 if disk edge is hit, otherwise 0.0:
        if (trackDiskHit == true)
        {
            image[index] = 0.0;
            pixelX[index] = xyz0[1];
            pixelY[index] = xyz0[2];
            pixelZ[index] = xyz0[3];
        }
    }
    
    // Save camera image:
    std::ofstream fileOut0(OUTPUTDIR + "Image" + metric.Name() + ".csv");
    fileOut0 << "#x,y,z,value\n";
    for(int ij=0; ij<camera.pixelCount; ij++)
        fileOut0 << pixelX[ij] << "," << pixelY[ij] << "," << pixelZ[ij] << "," << image[ij] << "\n";
    fileOut0.close();
}



int main()
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    int n = 10;
    int nx = n*21+1;
    int ny = n*27+1;
    int nz = n*12+1;
    Coord start(-14,-14,-1);
    Coord   end( 14, 22, 15);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    // Minkowski metric(grid, m, a);
    SchwarzSchild metricSS(grid, m, a);
    KerrSchild metricKerr(grid, m, a);
    // LebedevStencil lebedevStencil(5);

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

    // BlackHoleCsv(metric);
    // ThinDiskCsv(metric,diskInner,diskOuter);
    RaysFromCamera(metricSS,camera,diskInner,diskOuter, 1150);
    RaysFromCamera(metricKerr,camera,diskInner,diskOuter, 1200);
}