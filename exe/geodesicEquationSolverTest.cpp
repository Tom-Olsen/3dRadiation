#include <iostream>
#include "../src/GeodesicEquationSolver.h"
#include "../src/SpecialMath.h"
#include "../src/Stencil.h"
#include "../src/SphericalHarmonics.h"
using namespace std;


void StencilStreaming()
{
    ofstream fileOut((string)OUTPUTDIR + (string)"Test_StencilStreaming.txt");
    fileOut << "x, y, z, s \n";

    int n = 10; // grid size for test stencils
    size_t nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-4,-4,-4);
    Coord end(4,4,4);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);
    SchwarzSchild metric(grid, 1.0, 0.0);
    LebedevStencil stencil(21);
    double dt = 0.25 * (end[1] - start[1]) / n;
    cout << "Initialization complete." << endl;

    // Compute reversed geodesic for every direction in the stencil, at every test point:
    for(int k=0; k<n; k++)
    for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
    {
        double s = 1.0;
        Coord x0(start[1] + (i+0.5) * (end[1] - start[1]) / n, start[2] + (j+0.5) * (end[2] - start[2]) / n, start[3] + (k+0.5) * (end[3] - start[3]) / n);
        

        if (!metric.InsideBH(x0))
        {
            fileOut << x0[1] << ", " << x0[2] << ", " << x0[3] << ", " << 1.0 << "\n";
            float alpha = metric.GetAlpha(x0);

            for(int d=0; d<stencil.nDir; d++)
            {
                // Initial data for geodesic equation:
                double s = 1;
                Coord x = x0;
                Tensor3 c = stencil.Ct3(d);
                Tensor4 uIF(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
                Tensor3 vLF = Vec3ObservedByEulObs<IF, LF>(uIF, x, metric);

                // Solve geodesic equation backwards:
                s *= RK45_GeodesicEquation<-1>(dt, x, vLF, metric);
                fileOut << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
            }
        }
        cout << "Stencil(" << i << "," << j << "," << k << ") complete." << endl;
    }

    // Black Hole Horizon:
    double r = 2.0 * metric.m;
    for(int j=0; j<40; j++)
    for(int i=0; i<20; i++)
    {
        double theta = M_PI * (i + 0.5) / 20.0;
        double phi = 2.0 * M_PI * j / 39.0;
        double x = r * MySin(theta) * MyCos(phi);
        double y = r * MySin(theta) * MySin(phi);
        double z = r * MyCos(theta);
        fileOut << x << ", " << y << ", " << z << ", " << 1.0 << "\n";
    }


    fileOut.close();
}



void StencilStreamingWithHarmonics()
{
    ofstream fileOut0((string)OUTPUTDIR + (string)"Test_StencilStreamingWithHarmonicsBaseData.txt");
    ofstream fileOut1((string)OUTPUTDIR + (string)"Test_StencilStreamingWithHarmonics.txt");
    fileOut0 << "x, y, z, s \n";
    fileOut1 << "x, y, z, s \n";

    int n = 10; // grid size for test stencils
    size_t nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-4,-4,-4);
    Coord end(4,4,4);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);
    SchwarzSchild metric(grid, 1.0, 0.0);
    LebedevStencil stencil(21);
    LebedevStencil streamingStencil(5);
    double dt = 0.25 * (end[1] - start[1]) / n;
    cout << "Initialization complete." << endl;

    // Spherical Harmonics Coefficients:
    RealBuffer coefficientsS;   coefficientsS.resize(n * n * n * streamingStencil.nCoefficients);
    RealBuffer coefficientsX;   coefficientsX.resize(n * n * n * streamingStencil.nCoefficients);
    RealBuffer coefficientsY;   coefficientsY.resize(n * n * n * streamingStencil.nCoefficients);
    RealBuffer coefficientsZ;   coefficientsZ.resize(n * n * n * streamingStencil.nCoefficients);
    RealBuffer coefficientsCx;  coefficientsCx.resize(n * n * n * streamingStencil.nCoefficients);
    RealBuffer coefficientsCy;  coefficientsCy.resize(n * n * n * streamingStencil.nCoefficients);
    RealBuffer coefficientsCz;  coefficientsCz.resize(n * n * n * streamingStencil.nCoefficients);

    // Compute Spherical Harmonics Coefficients from streaming stencil:
    for(int k=0; k<n; k++)
    for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
    {
        double s = 1.0;
        Coord x0(start[1] + (i+0.5) * (end[1] - start[1]) / n, start[2] + (j+0.5) * (end[2] - start[2]) / n, start[3] + (k+0.5) * (end[3] - start[3]) / n);
        
        if (!metric.InsideBH(x0))
        {
            fileOut0 << x0[1] << ", " << x0[2] << ", " << x0[3] << ", " << 1.0 << "\n";
            float alpha = metric.GetAlpha(x0);

            double dataS[streamingStencil.nDir];
            double dataX[streamingStencil.nDir];
            double dataY[streamingStencil.nDir];
            double dataZ[streamingStencil.nDir];
            double dataCx[streamingStencil.nDir];
            double dataCy[streamingStencil.nDir];
            double dataCz[streamingStencil.nDir];
            
            for(int d=0; d<streamingStencil.nDir; d++)
            {
                // Initial data for geodesic equation:
                double s = 1;
                Coord x = x0;
                Tensor3 c = streamingStencil.Ct3(d);
                Tensor4 uIF(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
                Tensor3 vLF = Vec3ObservedByEulObs<IF, LF>(uIF, x, metric);
    
                // Solve geodesic equation backwards:
                s *= RK45_GeodesicEquation<-1>(dt, x, vLF, metric);
                fileOut0 << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
                dataS[d] = 1.0 / s;
                dataX[d] = x[1];
                dataY[d] = x[2];
                dataZ[d] = x[3];
                dataCx[d] = c[1];
                dataCy[d] = c[2];
                dataCz[d] = c[3];
                
                // Compute Spherical Harmonic Coefficinets:
                size_t ijk = i + j * n + k * n * n;
                size_t harmonicIndex = 0 + ijk * streamingStencil.nCoefficients;
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataS, &coefficientsS[harmonicIndex]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataX, &coefficientsX[harmonicIndex]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataY, &coefficientsY[harmonicIndex]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataZ, &coefficientsZ[harmonicIndex]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[harmonicIndex]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[harmonicIndex]);
                SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCz, &coefficientsCz[harmonicIndex]);
            }
        }
    }
    
    // Use spherical harmonics coefficients to compute origin point for every direction in the stencil:
    for(int k=0; k<n; k++)
    for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
    {
        size_t ijk = i + j * n + k * n * n;
        Coord x0(start[1] + (i+0.5) * (end[1] - start[1]) / n, start[2] + (j+0.5) * (end[2] - start[2]) / n, start[3] + (k+0.5) * (end[3] - start[3]) / n);
        if (!metric.InsideBH(x0))
        {
            fileOut1 << x0[1] << ", " << x0[2] << ", " << x0[3] << ", " << 1.0 << "\n";
            for(int d=0; d<stencil.nDir; d++)
            {
                size_t harmonicIndex = 0 + ijk * streamingStencil.nCoefficients;
                Coord x;
                Tensor3 c = stencil.Ct3(d);
                x[1] = SphericalHarmonicsXyz::GetValue(c, &coefficientsX[harmonicIndex], streamingStencil.nCoefficients);
                x[2] = SphericalHarmonicsXyz::GetValue(c, &coefficientsY[harmonicIndex], streamingStencil.nCoefficients);
                x[3] = SphericalHarmonicsXyz::GetValue(c, &coefficientsZ[harmonicIndex], streamingStencil.nCoefficients);
                double s = SphericalHarmonicsXyz::GetValue(c, &coefficientsS[harmonicIndex], streamingStencil.nCoefficients);
                fileOut1 << x[1] << ", " << x[2] << ", " << x[3] << ", " << s << "\n";
            }
        }
        cout << "Stencil(" << i << "," << j << "," << k << ") complete." << endl;
    }

    // Black Hole Horizon:
    double r = 2.0 * metric.m;
    for(int j=0; j<40; j++)
    for(int i=0; i<20; i++)
    {
        double theta = M_PI * (i + 0.5) / 20.0;
        double phi = 2.0 * M_PI * j / 39.0;
        double x = r * MySin(theta) * MyCos(phi);
        double y = r * MySin(theta) * MySin(phi);
        double z = r * MyCos(theta);
        fileOut0 << x << ", " << y << ", " << z << ", " << 1.0 << "\n";
        fileOut1 << x << ", " << y << ", " << z << ", " << 1.0 << "\n";
    }
    fileOut0.close();
    fileOut1.close();
}



void TestManyPhotons()
{
    ofstream fileOut((string)OUTPUTDIR + (string)"Test_GeodesicEquationSolver.txt");
    fileOut << "x, y, z, s \n";

    size_t nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-4,-4,-4);
    Coord end(4,4,4);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    cout << "Initialization complete." << endl;

    // Geodesics:
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

    // Black Hole Horizon:
    double r = 2.0 * metric.m;
    for(int j=0; j<20; j++)
    for(int i=0; i<10; i++)
    {
        double theta = M_PI * (i + 0.5) / 10.0;
        double phi = 2.0 * M_PI * j / 19.0;
        double x = r * MySin(theta) * MyCos(phi);
        double y = r * MySin(theta) * MySin(phi);
        double z = r * MyCos(theta);
        fileOut << x << ", " << y << ", " << z << ", " << 0.0 << "\n";
    }


    fileOut.close();
}



void BoundingGeodesics()
{
    ofstream fileOut0(OUTPUTDIR + (string)"ExactGeodesic0.txt");
    ofstream fileOut1(OUTPUTDIR + (string)"ExactGeodesic1.txt");
    fileOut0 << "x, y, z, s \n";
    fileOut1 << "x, y, z, s \n";

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
    // StencilStreaming();
    StencilStreamingWithHarmonics();
    // TestManyPhotons();
    // BoundingGeodesics();
}