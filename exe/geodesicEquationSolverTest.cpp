#include <iostream>
#include "../src/GeodesicEquationSolver.h"
#include "../src/SpecialMath.h"
using namespace std;


void TestManyPhotons()
{
    ofstream fileOut((string)OUTPUTDIR + (string)"Test_GeodesicEquationSolver.txt");
    fileOut << "#x, y, z, s \n";

    size_t nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-4,-4,-4);
    Coord end(4,4,4);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    cout << "Initialization complete." << endl;

    int n = 10;
    for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
    {
        double s = 1.0;
        Coord x(start[1] + (i+0.5) * (end[1] - start[1]) / n, start[2] + (j+0.5) * (end[2] - start[2]) / n, start[3]);
        Tensor4x4 g_ll = metric.GetMetric_ll(x);
        Tensor3x3 gamma_ll = metric.GetGamma_ll(x);
        Tensor3x3 delta_ll = metric.GetMinkowskiGamma_ll(x);

        Tensor4 uLF(1,0,0,1);
        uLF = NullNormalize(uLF, g_ll);
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF,x,metric);

        fileOut << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
        while(true)
        {
            if (grid.OutsideDomain(x) || metric.InsideBH(x))
                break;
            s *= RK45_GeodesicEquation<1>(5 * grid.dt, x, vLF, metric);
            fileOut << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
        }
        cout << "Photon(" << i << "," << j << ") complete." << endl;
    }

    fileOut.close();
}



void BoundingGeodesics()
{
    ofstream fileOut0(OUTPUTDIR + (string)"ExactGeodesic0.txt");
    ofstream fileOut1(OUTPUTDIR + (string)"ExactGeodesic1.txt");
    fileOut0 << "#x, y, z, s \n";
    fileOut1 << "#x, y, z, s \n";

    size_t nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-0.0001,0,-0.5);
    Coord end(5,4,0.5);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    cout << "Initialization complete." << endl;

    // First Photon:
    {
        double s = 1.0;
        Coord x(0.0, 3.0, 0.0);
        Tensor4x4 g_ll = metric.GetMetric_ll(x);
        Tensor3x3 gamma_ll = metric.GetGamma_ll(x);
        Tensor3x3 delta_ll = metric.GetMinkowskiGamma_ll(x);

        Tensor4 uLF(1,1,0,0);
        uLF = NullNormalize(uLF, g_ll);
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF,x,metric);

        fileOut0 << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
        while(true)
        {
            if (grid.OutsideDomain(x) || metric.InsideBH(x))
                break;
            s *= RK45_GeodesicEquation<1>(5 * grid.dt, x, vLF, metric);
            fileOut0 << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
        }
        cout << "Photon1 complete." << endl;
        fileOut0.close();
    }
    // Second Photon:
    {
        double s = 1.0;
        Coord x(0.0, 3.5, 0.0);
        Tensor4x4 g_ll = metric.GetMetric_ll(x);
        Tensor3x3 gamma_ll = metric.GetGamma_ll(x);
        Tensor3x3 delta_ll = metric.GetMinkowskiGamma_ll(x);

        Tensor4 uLF(1,1,0,0);
        uLF = NullNormalize(uLF, g_ll);
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF,x,metric);

        fileOut1 << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
        while(true)
        {
            if (grid.OutsideDomain(x) || metric.InsideBH(x))
                break;
            s *= RK45_GeodesicEquation<1>(5 * grid.dt, x, vLF, metric);
            fileOut1 << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
        }
        cout << "Photon2 complete." << endl;
        fileOut1.close();
    }

}



int main()
{
    // TestManyPhotons();
    BoundingGeodesics();
}