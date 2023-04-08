#include <iostream>
#include "../src/GeodesicEquationSolver.h"
#include "../src/TensorOperations.h"
using namespace std;



int main()
{
    std::ofstream fileOut((std::string)OUTPUTDIR + (std::string)"Test_GeodesicEquationSolver.csv");
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
        Tensor3x3 delta_ll(1,0,0, 0,1,0, 0,0,1);

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