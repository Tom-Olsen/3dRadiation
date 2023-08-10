#include <iostream>
#include "../src/Radiation.h"
using namespace std;



void SphericalHarmonicsExpansion()
{
    size_t nx, ny, nz;
    nx = ny = nz = 20;
    Coord start(0,0,0);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    GaussLegendreStencil stencil(15);
    LebedevStencil lebedevStencil(5);
    InterpolationGrid interpGrid(500, 1000, stencil);
    Camera camera;
    Config config =
    {
        .name = "bigFunctionsTest",
        .simTime = 1,
        .writePeriod = 1.0/100,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = false,
        .printToTerminal = false,
        .useCamera = false,
        .streamingType = StreamingType::CurvedAdaptive,
    };

    Radiation radiation(metric, stencil, lebedevStencil, interpGrid, camera, config);
    radiation.UpdateSphericalHarmonicsCoefficients();
    
    ofstream file0(OUTPUTDIR + "SphericalHarmonicsExpansionCoord.txt");
    ofstream file1(OUTPUTDIR + "SphericalHarmonicsExpansionVeloc.txt");
    ofstream file2(OUTPUTDIR + "GeodesicCoord.txt");
    ofstream file3(OUTPUTDIR + "GeodesicVeloc.txt");
    file0 << "#x, y, z, s \n";
    file1 << "#x, y, z, s \n";
    file2 << "#x, y, z, s \n";
    file3 << "#x, y, z, s \n";
	for(size_t k=2; k<grid.nz-2; k+=2)
	for(size_t j=2; j<grid.ny-2; j+=2)
	for(size_t i=2; i<grid.nx-2; i+=2)
    {
        size_t ijk = grid.Index(i,j,k);
        double alpha = metric.GetAlpha(ijk);
        Coord center = grid.xyz(i,j,k);
        for(size_t d=0; d<stencil.nDir; d++)
        {
            // Get pos and vel by Spherical Harmonic Expansion:
            {
                Tensor3 direction = stencil.Ct3(d);
                double s = radiation.GetFrequencyShift(ijk,direction);
                Coord xyz = radiation.GetTempCoordinate(ijk,direction);
                Tensor3 vIF = radiation.GetTemp3VelocityIF(ijk,direction);
                if(!metric.InsideBH(xyz))
                {
                    file0 << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << s << "\n";
                    file1 << center[1] + 0.1*vIF[1] << ", " << center[2] + 0.1*vIF[2] << ", " << center[3] + 0.1*vIF[3] << ", " << s << "\n";
                }
            }
            
            // Get pos and vel by Geodesic Equation Solver:
            {
                double s = 1;
			    Coord xyz = center;
                Tensor3 c = stencil.Ct3(d);
                Tensor4 u(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
                Tensor3 vLF = Vec3ObservedByEulObs<IF,LF>(u, xyz, metric);
    
			    if(!metric.InsideBH(xyz))
                {
    	        	s *= RK45_GeodesicEquation<-1>(grid.dt, xyz, vLF, metric);
                    Tensor3 vIF = TransformLFtoIF(vLF, metric.GetTetradInverse(xyz));
                    file2 << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << 1.0/s << "\n";
                    file3 << center[1] + 0.1*vIF[1] << ", " << center[2] + 0.1*vIF[2] << ", " << center[3] + 0.1*vIF[3] << ", " << 1.0/s << "\n";
                }
            }
        }
    }
    file0.close();
    file1.close();
    file2.close();
    file3.close();

    cout << "4 files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



int main()
{
    SphericalHarmonicsExpansion();
}