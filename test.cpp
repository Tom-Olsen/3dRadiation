#include <iostream>
#include "src/Includes.hh"

using namespace std;


void Test_TensorTypes()
{
    // Coord:
    Coord xyz(1.1,2.2,3.3);
    xyz.Print("xyz");
    PrintDouble(xyz[1],"x");
    PrintDouble(xyz[2],"y");
    PrintDouble(xyz[3],"z");
    cout << endl;

    // Tensor3:
    cout << "Tensor3" << endl;
    Tensor3 a3(1,2,3);
    a3.Print("a3");
    for(int i=1; i<4; i++)
        PrintDouble(a3[i],"a3[" + std::to_string(i) + "]");
    Tensor3 b3(4); 
    b3.Print("b3");
    for(int i=1; i<4; i++)
        PrintDouble(b3[i],"a3[" + std::to_string(i) + "]");
    cout << endl;
    
    // Tensor4:
    cout << "Tensor4" << endl;
    Tensor4 a4(1,2,3,4);
    a4.Print("a4");
    for(int i=0; i<4; i++)
        PrintDouble(a4[i],"a4[" + std::to_string(i) + "]");
    Tensor4 b4(4); 
    b4.Print("b4");
    for(int i=0; i<4; i++)
        PrintDouble(b4[i],"b4[" + std::to_string(i) + "]");
    cout << endl;

    // Tensor3x3:
    cout << "Tensor3x3" << endl;
    Tensor3x3 A33(0,1,2, 3,4,5, 6,7,8);
    A33.Print("A33");
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            PrintDouble(A33[{i,j}],"A33[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor3x3 B33(9);
    B33.Print("A33");
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            PrintDouble(B33[{i,j}],"A33[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;

    // Tensor4x4:
    cout << "Tensor4x4" << endl;
    Tensor4x4 A44(0,1,2,3, 2,3,4,5, 3,4,5,6, 4,5,6,7);
    A44.Print("A44");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            PrintDouble(A44[{i,j}],"A44[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor4x4 B44(1.234);
    B44.Print("A44");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            PrintDouble(B44[{i,j}],"A44[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;

    // Tensor3x3x3:
    cout << "Tensor3x3x3" << endl;
    Tensor3x3x3 A333(0);
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            for(int k=1; k<4; k++)
                A333[{i,j,k}] = (k-1) + 3*(j-1) + 9*(i-1);
    A333.Print("A333");
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            for(int k=1; k<4; k++)
                PrintDouble(A333[{i,j,k}],"A333[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;

    // Tensor4x4x4:
    cout << "Tensor4x4x4" << endl;
    Tensor4x4x4 A444(0);
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                A444[{i,j,k}] = k + 4*j + 16*i;
    A444.Print("A444");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                PrintDouble(A444[{i,j,k}],"A444[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;
}



void Test_Grid()
{
    size_t nx, ny, nz;
    nx = ny = nz = 20;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    auto dist = [] (Coord a, Coord b)
    {
        double x = (a[1] - b[1]);
        double y = (a[2] - b[2]);
        double z = (a[3] - b[3]);
        return sqrt(x*x + y*y + z*z);
    };

    Coord rOrb(0.3,0.3,0.3);
    Coord gOrb(0.3,0.7,0.4);
    Coord bOrb(0.7,0.4,0.6);
    Coord aOrb(0.6,0.7,0.6);
    RealBuffer data0(grid.nxyz);
    RealBuffer data1(grid.nxyz);
    RealBuffer data2(grid.nxyz);
    RealBuffer data3(grid.nxyz);
    for(size_t i=0; i<grid.nx; i++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t k=0; k<grid.nz; k++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        data0[ijk] = (dist(rOrb,xyz) < 0.2) ? 1 : 0;
        data1[ijk] = (dist(gOrb,xyz) < 0.2) ? 1 : 0;
        data2[ijk] = (dist(bOrb,xyz) < 0.2) ? 1 : 0;
        data3[ijk] = 1;
    }
    // grid.WriteFrametoJson(0,data0,data1,data2,data3,0,"output","Test_Grid");
    grid.WriteFrametoCsv (0,data0,data1,data2,data3,0,"output","Test_Grid");
}



void Test_Interpolation()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    RealBuffer data0(grid.nxyz);
    RealBuffer data1(grid.nxyz);
    RealBuffer data2(grid.nxyz);
    RealBuffer data3(grid.nxyz);
    for(size_t i=0; i<grid.nx; i++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t k=0; k<grid.nz; k++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        data0[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 1,0,0,0,0,0,0,1);
        data1[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,1,0,0,0,0,1,0);
        data2[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,0,1,0,0,1,0,0);
        data3[ijk] = 1;
    }
    grid.WriteFrametoCsv(0,data0,data1,data2,data3,0,"output","Test_Interpolation");
}



void Test_Metric()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(1,1,1);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);

    Coord xyz(1.5, 2, 2.5);
    Tensor4x4 g_ll = metric.GetMetric_ll(xyz);
    Tensor4x4 tetrad = metric.GetTetrad(xyz);
    Tensor4x4 eta(0);

    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
    for(int I=0; I<4; I++)
    for(int J=0; J<4; J++)
        eta[{i,j}] += tetrad[{I,i}] * tetrad[{J,j}] * g_ll[{I,J}];

    g_ll.Print("  g_ll",1);
    tetrad.Print("tetrad",1);
    eta.Print("   eta",1);
}



void Test_GeodesicEquationSolver()
{
    std::ofstream fileOut("output/Test_GeodesicEquationSolver.csv");
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



void Test_Stencil()
{
    Stencil stencil(5);
    stencil.Print();

    {
        double d0 = 2.3;
        double d1 = 6.7;
        double x = stencil.Cx(d0,d1);
        double y = stencil.Cy(d0,d1);
        double z = stencil.Cz(d0,d1);
        double theta = stencil.Theta(d0,d1);
        double phi = stencil.Phi(d0,d1);
        Tensor3 cxyz = stencil.Cxyz(d0,d1);
        cout << theta << ", " << cxyz.Theta() << endl;
        cout << phi << ", " << cxyz.Phi() << endl;
    }
    cout << endl;
}
void Test_LebedevStencil()
{
    LebedevStencil3 lebedev3;
    lebedev3.Print();
    LebedevStencil5 lebedev5;
    lebedev5.Print();
    LebedevStencil7 lebedev7;
    lebedev7.Print();
}



void Test_Quadrature()
{
    size_t nTh = 19;
    Stencil stencil(nTh);

    double data[stencil.nThPh];
    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        double x = stencil.Cx(d0,d1);
        double y = stencil.Cy(d0,d1);
        double z = stencil.Cz(d0,d1);

        size_t d = stencil.Index(d0,d1);
        data[d] = 0.1 *SphericalHarmonics::Y00(x,y,z);

        data[d] += 1.0 * SphericalHarmonics::Y1m1(x,y,z);
        data[d] += 1.1 * SphericalHarmonics::Y10(x,y,z);
        data[d] += 1.2 * SphericalHarmonics::Y1p1(x,y,z);
        
        data[d] += 2.1 * SphericalHarmonics::Y2m2(x,y,z);
        data[d] += 2.2 * SphericalHarmonics::Y2m1(x,y,z);
        data[d] += 2.3 * SphericalHarmonics::Y20(x,y,z);
        data[d] += 2.4 * SphericalHarmonics::Y2p1(x,y,z);
        data[d] += 2.5 * SphericalHarmonics::Y2p2(x,y,z);
        
        data[d] += 3.0 * SphericalHarmonics::Y3m3(x,y,z);
        data[d] += 3.1 * SphericalHarmonics::Y3m2(x,y,z);
        data[d] += 3.2 * SphericalHarmonics::Y3m1(x,y,z);
        data[d] += 3.3 * SphericalHarmonics::Y30 (x,y,z);
        data[d] += 3.4 * SphericalHarmonics::Y3p1(x,y,z);
        data[d] += 3.5 * SphericalHarmonics::Y3p2(x,y,z);
        data[d] += 3.6 * SphericalHarmonics::Y3p3(x,y,z);
        
        data[d] += 4.0 * SphericalHarmonics::Y4m4(x,y,z);
        data[d] += 4.1 * SphericalHarmonics::Y4m3(x,y,z);
        data[d] += 4.2 * SphericalHarmonics::Y4m2(x,y,z);
        data[d] += 4.3 * SphericalHarmonics::Y4m1(x,y,z);
        data[d] += 4.4 * SphericalHarmonics::Y40 (x,y,z);
        data[d] += 4.5 * SphericalHarmonics::Y4p1(x,y,z);
        data[d] += 4.6 * SphericalHarmonics::Y4p2(x,y,z);
        data[d] += 4.7 * SphericalHarmonics::Y4p3(x,y,z);
        data[d] += 4.8 * SphericalHarmonics::Y4p4(x,y,z);
    }
    vector<double>c = SphericalHarmonics::GetCoefficients(stencil,data,min(stencil.nCoefficients,(size_t)25));
    
    for(size_t i=0; i<c.size(); i++)
        cout << "c" << i << ": " << c[i] << endl;
    cout << endl;

    cout << "Real Value, Compressed Value:" << endl;
    size_t d0 = stencil.nTh/2;
    size_t d1 = stencil.nPh/2;
    double theta = stencil.Theta(d0,d1);
    double phi = stencil.Phi(d0,d1);
    cout << data[stencil.Index(d0,d1)] << ", " << SphericalHarmonics::GetValue(theta,phi,c,min(stencil.nCoefficients,(size_t)25)) << endl << endl;

    cout << "nCoefficients = " << stencil.nCoefficients << endl;
    // stencil.Print();
    stencil.WriteToCsv();
}
void Test_QuadratureLebedev()
{
    LebedevStencil7 lebedev;

    double data[lebedev.nDir];
    for(size_t d=0; d<lebedev.nDir; d++)
    {
        double x = lebedev.Cx(d);
        double y = lebedev.Cy(d);
        double z = lebedev.Cz(d);

        data[d] = 0.1 *SphericalHarmonics::Y00(x,y,z);

        data[d] += 1.0 * SphericalHarmonics::Y1m1(x,y,z);
        data[d] += 1.1 * SphericalHarmonics::Y10(x,y,z);
        data[d] += 1.2 * SphericalHarmonics::Y1p1(x,y,z);
        
        data[d] += 2.1 * SphericalHarmonics::Y2m2(x,y,z);
        data[d] += 2.2 * SphericalHarmonics::Y2m1(x,y,z);
        data[d] += 2.3 * SphericalHarmonics::Y20(x,y,z);
        data[d] += 2.4 * SphericalHarmonics::Y2p1(x,y,z);
        data[d] += 2.5 * SphericalHarmonics::Y2p2(x,y,z);
        
        data[d] += 3.0 * SphericalHarmonics::Y3m3(x,y,z);
        data[d] += 3.1 * SphericalHarmonics::Y3m2(x,y,z);
        data[d] += 3.2 * SphericalHarmonics::Y3m1(x,y,z);
        data[d] += 3.3 * SphericalHarmonics::Y30 (x,y,z);
        data[d] += 3.4 * SphericalHarmonics::Y3p1(x,y,z);
        data[d] += 3.5 * SphericalHarmonics::Y3p2(x,y,z);
        data[d] += 3.6 * SphericalHarmonics::Y3p3(x,y,z);
        
        data[d] += 4.0 * SphericalHarmonics::Y4m4(x,y,z);
        data[d] += 4.1 * SphericalHarmonics::Y4m3(x,y,z);
        data[d] += 4.2 * SphericalHarmonics::Y4m2(x,y,z);
        data[d] += 4.3 * SphericalHarmonics::Y4m1(x,y,z);
        data[d] += 4.4 * SphericalHarmonics::Y40 (x,y,z);
        data[d] += 4.5 * SphericalHarmonics::Y4p1(x,y,z);
        data[d] += 4.6 * SphericalHarmonics::Y4p2(x,y,z);
        data[d] += 4.7 * SphericalHarmonics::Y4p3(x,y,z);
        data[d] += 4.8 * SphericalHarmonics::Y4p4(x,y,z);
    }
    vector<double>c = SphericalHarmonics::GetCoefficients(lebedev,data,min(lebedev.nCoefficients,(size_t)25));

    for(size_t i=0; i<c.size(); i++)
        cout << "c" << i << ": " << c[i] << endl;
    cout << endl;
    
    cout << "Real Value, Compressed Value:" << endl;
    size_t d = lebedev.nDir / 2;
    double theta = lebedev.Theta(d);
    double phi = lebedev.Phi(d);
    cout << data[d] << ", " << SphericalHarmonics::GetValue(theta,phi,c,min(lebedev.nCoefficients,(size_t)25)) << endl << endl;

    cout << "         nDir = " << lebedev.nDir << endl;
    cout << "       nOrder = " << lebedev.nOrder << endl;
    cout << "nCoefficients = " << lebedev.nCoefficients << endl;
    
    lebedev.Print();
    lebedev.WriteToCsv();
}
void Test_QuadratureGaussLegendre()
{
    // SphericalHarmonics::GetCoefficients only supported up to 9
    GaussLegendreStencil19 gauss;

    double data[gauss.nDir];
    for(size_t d=0; d<gauss.nDir; d++)
    {
        double x = gauss.Cx(d);
        double y = gauss.Cy(d);
        double z = gauss.Cz(d);

        data[d] = 0.1 *SphericalHarmonics::Y00(x,y,z);

        data[d] += 1.0 * SphericalHarmonics::Y1m1(x,y,z);
        data[d] += 1.1 * SphericalHarmonics::Y10(x,y,z);
        data[d] += 1.2 * SphericalHarmonics::Y1p1(x,y,z);
        
        data[d] += 2.1 * SphericalHarmonics::Y2m2(x,y,z);
        data[d] += 2.2 * SphericalHarmonics::Y2m1(x,y,z);
        data[d] += 2.3 * SphericalHarmonics::Y20(x,y,z);
        data[d] += 2.4 * SphericalHarmonics::Y2p1(x,y,z);
        data[d] += 2.5 * SphericalHarmonics::Y2p2(x,y,z);
        
        data[d] += 3.0 * SphericalHarmonics::Y3m3(x,y,z);
        data[d] += 3.1 * SphericalHarmonics::Y3m2(x,y,z);
        data[d] += 3.2 * SphericalHarmonics::Y3m1(x,y,z);
        data[d] += 3.3 * SphericalHarmonics::Y30 (x,y,z);
        data[d] += 3.4 * SphericalHarmonics::Y3p1(x,y,z);
        data[d] += 3.5 * SphericalHarmonics::Y3p2(x,y,z);
        data[d] += 3.6 * SphericalHarmonics::Y3p3(x,y,z);
        
        data[d] += 4.0 * SphericalHarmonics::Y4m4(x,y,z);
        data[d] += 4.1 * SphericalHarmonics::Y4m3(x,y,z);
        data[d] += 4.2 * SphericalHarmonics::Y4m2(x,y,z);
        data[d] += 4.3 * SphericalHarmonics::Y4m1(x,y,z);
        data[d] += 4.4 * SphericalHarmonics::Y40 (x,y,z);
        data[d] += 4.5 * SphericalHarmonics::Y4p1(x,y,z);
        data[d] += 4.6 * SphericalHarmonics::Y4p2(x,y,z);
        data[d] += 4.7 * SphericalHarmonics::Y4p3(x,y,z);
        data[d] += 4.8 * SphericalHarmonics::Y4p4(x,y,z);
    }
    vector<double>c = SphericalHarmonics::GetCoefficients(gauss,data,min(gauss.nCoefficients,(size_t)25));

    for(size_t i=0; i<c.size(); i++)
        cout << "c" << i << ": " << c[i] << endl;
    cout << endl;
    
    cout << "Real Value, Compressed Value:" << endl;
    size_t d = gauss.nDir / 2;
    double theta = gauss.Theta(d);
    double phi = gauss.Phi(d);
    cout << data[d] << ", " << SphericalHarmonics::GetValue(theta,phi,c,min(gauss.nCoefficients,(size_t)25)) << endl << endl;

    cout << "         nDir = " << gauss.nDir << endl;
    cout << "       nOrder = " << gauss.nOrder << endl;
    cout << "nCoefficients = " << gauss.nCoefficients << endl;
    
    // gauss.Print();
    gauss.WriteToCsv();
}



void Test_SphericalHarmonicsExpansion()
{
    size_t nx, ny, nz;
    nx = ny = nz = 20;
    Coord start(0,0,0);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    Stencil stencil(15);
    LebedevStencil5 lebedevStencil;
    Camera camera;

    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedDynamic);
    radiation.UpdateSphericalHarmonicsCoefficients();
    
    ofstream file0("output/Test_SphericalHarmonicsExpansionCoord.csv");
    ofstream file1("output/Test_SphericalHarmonicsExpansionVeloc.csv");
    ofstream file2("output/Test_GeodesicCoord.csv");
    ofstream file3("output/Test_GeodesicVeloc.csv");
    file0 << "#x, y, z, s \n";
    file1 << "#x, y, z, s \n";
    file2 << "#x, y, z, s \n";
    file3 << "#x, y, z, s \n";
	for(size_t k=2; k<grid.nz-2; k+=2)
	for(size_t j=2; j<grid.ny-2; j+=2)
	for(size_t i=2; i<grid.nx-2; i+=2)
    {
		if(i == grid.nx/2 && j == grid.ny/2 && k == grid.nz/2)
        {
	        std::vector<double> cS(lebedevStencil.nCoefficients);
	        for(size_t f=0; f<lebedevStencil.nCoefficients; f++)
	        	cS[f] = radiation.coefficientsS[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
                
    		for(size_t i=0; i<cS.size(); i++)
    		    std::cout << "cS" << i << ": " << cS[i] << std::endl;
    		std::cout << std::endl;
        }

        size_t ijk = grid.Index(i,j,k);
        double alpha = metric.GetAlpha(ijk);
        Coord center = grid.xyz(i,j,k);
        for(size_t d1=0; d1<stencil.nPh; d1++)
        for(size_t d0=0; d0<stencil.nTh; d0++)
        {
            // Get pos and vel by Spherical Harmonic Expansion:
            {
                Tensor3 direction = stencil.Cxyz(d0,d1);
                double s = radiation.GetFrequencyShift(ijk,direction);
                Coord xyz = radiation.GetTempCoordinate(ijk,direction);
                Tensor3 v = radiation.GetTemp3Velocity(ijk,direction);
                if(!metric.InsideBH(xyz))
                {
                    file0 << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << s << "\n";
                    file1 << center[1] + 0.1*v[1] << ", " << center[2] + 0.1*v[2] << ", " << center[3] + 0.1*v[3] << ", " << s << "\n";
                }
            }
            
            // Get pos and vel by Geodesic Equation Solver:
            {
                double s = 1;
			    Coord xyz = center;
                Tensor3 c = stencil.Cxyz(d0,d1);
                Tensor4 u(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
                Tensor3 v = Vec3ObservedByEulObs<IF,LF>(u, xyz, metric);
    
			    if(!metric.InsideBH(xyz))
                {
    	        	s *= RK45_GeodesicEquation<-1>(grid.dt, xyz, v, metric);
                    file2 << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << 1.0/s << "\n";
                    file3 << center[1] + 0.1*v[1] << ", " << center[2] + 0.1*v[2] << ", " << center[3] + 0.1*v[3] << ", " << 1.0/s << "\n";
                }
            }
        }
    }
    file0.close();
    file1.close();
    file2.close();
    file3.close();
}



void Test_QuaternionRotation()
{
    using namespace glm;
    vec3 from(0,0,1);
    vec3 to(1,0,0);
    quat q = quat(from, to);
    // quat qInv = glm::conjugate(q);
    quat qInv = quat(to, from);

    Stencil stencil(15);

    ofstream file0("output/Test_QuaternionRotation Normal");
    ofstream file1("output/Test_QuaternionRotation Rotated");
    file0 << "#x, y, z" << endl;
    file1 << "#x, y, z" << endl;
    
    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        Tensor3 p = stencil.Cxyz(d0,d1);
        file0 << p[1] << ", " << p[2] << ", " << p[3] << endl;
        p = q * p;
        file1 << p[1] << ", " << p[2] << ", " << p[3] << endl;
    }

    file0.close();
    file1.close();

    Tensor3 p(1,2,3);
    Tensor3 pRot = q * p;
    Tensor3 pRotInvRot = qInv * pRot;

    p.Print("   p");
    pRot.Print("pRot");
    pRotInvRot.Print("  p?");
}



void Test_IntensityAt()
{
    // This test is simulating the Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF) function.
    Stencil stencil(5);
    LebedevStencil11 lebedev;

    double I[stencil.nThPh];
    double Inorth;
    double Isouth;

    ofstream file0("output/Test_IntensityAtStencil");
    ofstream file1("output/Test_IntensityAtLebedev");
    file0 << "x, y, z, color\n";
    file1 << "x, y, z, color\n";

    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        size_t d = stencil.Index(d0,d1);
        Tensor3 xyz = stencil.Cxyz(d0,d1);
        I[d] = 1.0 + 0.75 * SphericalHarmonics::Y1p1(xyz) + 0.75 * SphericalHarmonics::Y10 (xyz) + 0.75 * SphericalHarmonics::Y1m1(xyz);
        file0 << I[d] * xyz[1] << ", " << I[d] * xyz[2] << ", " << I[d] * xyz[3] << ", " << 1 << "\n";
    }
    Inorth = 1.0 + 0.75 * SphericalHarmonics::Y1p1(0, 0, 1) + 0.75 * SphericalHarmonics::Y10 (0, 0, 1) + 0.75 * SphericalHarmonics::Y1m1(0, 0, 1);
    Isouth = 1.0 + 0.75 * SphericalHarmonics::Y1p1(0, 0,-1) + 0.75 * SphericalHarmonics::Y10 (0, 0,-1) + 0.75 * SphericalHarmonics::Y1m1(0, 0,-1);
    file0 << Inorth * 0 << ", " << Inorth * 0 << ", " << Inorth * 1 << ", "  << 0 << "\n";
    file0 << Isouth * 0 << ", " << Isouth * 0 << ", " << Isouth * -1 << ", " << 0 << "\n";

    for(size_t d=0; d<lebedev.nDir; d++)
    {
        Tensor3 xyz = lebedev.Cxyz(d);
	    double t = stencil.d0(xyz.Theta());
	    double p = stencil.d1(xyz.Phi());
        
        // Don't use size_t due to possible negative values.
        int t0 = floor(t);
        int t1 = t0 + 1;
        int p0 = ((int)floor(p) + stencil.nPh) % stencil.nPh;
        int p1 = (p0 + 1) % stencil.nPh;
        
	    // Get intensities at nearest directions. Use north pole for t0==-1 and south pole for t1==nTh.
        double I00 = (t0 == -1)          ? Inorth : I[stencil.Index(t0,p0)];
        double I01 = (t0 == -1)          ? Inorth : I[stencil.Index(t0,p1)];
        double I10 = (t1 == stencil.nTh) ? Isouth : I[stencil.Index(t1,p0)];
        double I11 = (t1 == stencil.nTh) ? Isouth : I[stencil.Index(t1,p1)];

	    // Fractional part of sphere grid index. Near poles it must be streched:
	    // South: from [0.0, 0.5] to [0.0, 1.0]
	    // North: from [0.5, 1.0] to [0.0, 1.0]
	    double tfrac = t - floor(t);	// Bug Fixing: Seems to be always 0 or 1
	    double pfrac = p - floor(p);
        if (t0 == -1)			// North
            tfrac = 2.0 * (tfrac - 0.5);
        if (t1 == stencil.nTh)	// South
            tfrac *= 2.0;

	    double I = BilinearInterpolation(tfrac, pfrac, I00, I01, I10, I11);
        double Iexpected = 1.0 + 0.75 * SphericalHarmonics::Y1p1(xyz) + 0.75 * SphericalHarmonics::Y10 (xyz) + 0.75 * SphericalHarmonics::Y1m1(xyz);
        double color = 0;
        if(abs(I - Iexpected) > 0.1)
        {   
            PrintDouble(t    ,"    t");
            PrintDouble(t0   ,"   t0");
            PrintDouble(t1   ,"   t1");
            PrintDouble(tfrac,"tfrac",1);
            PrintDouble(p    ,"    p");
            PrintDouble(p0   ,"   p0");
            PrintDouble(p1   ,"   p1");
            PrintDouble(pfrac,"pfrac",1);
            color = 1;
        }
        file1 << I * xyz[1] << ", " << I * xyz[2] << ", " << I * xyz[3] << ", " << color << "\n";
    }
    file0.close();
    file1.close();
}



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
    Stencil stencil(nTh);
    stencil.sigma = sigma;
    LebedevStencil5 lebedevStencil;

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



void Test_MyAtan2()
{
    int n = 5;
    std::cout << "atan2f(y,x):" << std::endl;
    for(int j=0; j<n; j++)
    {
        float y = -2.0 * j / (n-1.0) + 1.0;
        for(int i=0; i<n; i++)
        {
            float x = 2.0 * i / (n-1.0) - 1.0;
            float atan2 = atan2f(y,x);
            std::cout << "(" << Format(x,3) << "," << Format(y,3) << "," << Format(atan2,3) << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    std::cout << "MyAtan2(y,x):" << std::endl;
    for(int j=0; j<n; j++)
    {
        float y = -2.0 * j / (n-1.0) + 1.0;
        for(int i=0; i<n; i++)
        {
            float x = 2.0 * i / (n-1.0) - 1.0;
            float atan2 = MyAtan2(y,x);
            std::cout << "(" << Format(x,3) << "," << Format(y,3) << "," << Format(atan2,3) << ") ";
        }
        std::cout << std::endl;
    }
}



void Test_MySin_MyCos()
{
    int n = 20;
    for(int i=0; i<n; i++)
    {
        float phi = 2 * M_PI * i / (n - 1.0);
        Tensor3(phi, sin(phi), MySin<9>(phi)).Print("(phi, sin, MySin)");
    }

    std::cout << std::endl << std::endl;
    
    for(int i=0; i<n; i++)
    {
        float phi = 2 * M_PI * i / (n - 1.0);
        Tensor3(phi, cos(phi), MyCos<9>(phi)).Print("(phi, cos, MyCos)");
    }
}



void Test_PhotonFourVelocity()
{
    // Grid, Metric, Stencil:
    Coord start(2,2,2);
    Coord end(3,3,3);
    Grid grid(50, 50, 50, start, end);
    SchwarzSchild metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);

    Coord x(2.2,2.3,2.4);
    double alpha = metric.GetAlpha(x);
    Tensor4x4 g_ll = metric.GetMetric_ll(x);

    Tensor4 u(1,0,0,1);
    u = NullNormalize(u,g_ll);
    double uNorm2 = Norm2(u,g_ll);
    Tensor4 h(1,0,0,alpha);
    double hNorm2 = Norm2(h,g_ll);

    u.Print("u");
    PrintDouble(uNorm2,"|u|^2");
    h.Print("h");
    PrintDouble(hNorm2,"|h|^2");
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
        double radius = xyz.Radius();

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
    Stencil stencil(7);
    stencil.sigma = 1.0;
    LebedevStencil3 lebedevStencil;
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);

    // Initial Data:
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double r = xyz.Radius();
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
    Stencil stencil(7);
    stencil.sigma = 1.0;
    LebedevStencil3 lebedevStencil;
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);

    // Initial Data:
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double r = xyz.Radius();
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
    Stencil stencil(15);
    stencil.sigma = 100;
    LebedevStencil3 lebedevStencil;
    Camera camera;
    // Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);

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



void StreamCurvedBeam(size_t nx, size_t ny, size_t nz, size_t nTh, int sigma, int simTime)
{
    // Grid, Metric, Stencil:
    Coord start(-1,0,-0.1);
    Coord end(1,3.6,3.9);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);   // needs at least LebedevStencil5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    Stencil stencil(nTh);
    stencil.sigma = sigma;
    LebedevStencil5 lebedevStencil;

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
            Tensor4 u(1,0,0,1);
            u = NullNormalize(u,metric.GetMetric_ll(ijk));
            Tensor3 v = Vec3ObservedByEulObs<LF,IF>(u, xyz, metric);

            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
            radiation.initialNx[ijk] = v[1];
            radiation.initialNy[ijk] = v[2];
            radiation.initialNz[ijk] = v[3];
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
        }
    }

    // Start simulation:
    Config config =
    {
        // .name = "Curved Beam Static " + metric.Name() + " " + std::to_string(stencil.nTh) + "." + std::to_string(stencil.nPh)
            //   + " s" + std::to_string(sigma) + " Leb" + std::to_string(lebedevStencil.nOrder) + " t" + std::to_string(simTime),
        .name = "Curved Beam Dynamic " + metric.Name() + " " + std::to_string(stencil.nTh) + "." + std::to_string(stencil.nPh)
              + " " + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z"
              + " s" + std::to_string(sigma) + " Leb" + std::to_string(lebedevStencil.nOrder) + " t" + std::to_string(simTime),
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



void ThinDisk(size_t nx, size_t ny, size_t nz, size_t nTh, int sigma, int simTime, bool diskIsHomogeneous)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    Coord start(-14,-14,-8);
    Coord   end( 14, 22, 12);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, m, a);
    Stencil stencil(nTh);
    stencil.sigma = sigma;
    LebedevStencil5 lebedevStencil;

    // Camera:
    size_t resX = 400;
    size_t resY = 200;
    size_t width = 26;
    size_t height = 13;
    Coord position(0,19,3);
    double degreeToRadians = 2.0 * M_PI / 360.0;
    double angleX = 100 * degreeToRadians;
    double angleY = 0 * degreeToRadians;
    double angleZ = 0 * degreeToRadians;
    glm::vec3 eulerAngles(angleX,angleY,angleZ);
    Camera camera(resX, resY, width, height, position, eulerAngles);

    // Radiation:
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedDynamic);

    // Initial Data:
    #pragma omp parallel for
    for(size_t k=0; k<grid.nz; k++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double radius = xyz.Radius();

        // Disk:
        if(diskInner <= radius && radius <= diskOuter && abs(xyz[3]) < 0.9 * grid.dz)
        {
            radiation.isInitialGridPoint[ijk] = false;
            radiation.initialE[ijk] = 0;
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            if (diskIsHomogeneous)
                radiation.initialEta[ijk] = 1;
            EIGEN_DEFAULT_SETTINGS_H
            {
                // s goes from 0 to 1.
                double d = (diskOuter - xyz[2]) / (2.0 * diskOuter);
                radiation.initialEta[ijk] = 0.1 + 0.9*d*d*d*d;
            }
        }
    }

    // Start simulation:
    Config config =
    {
        .name = "Thin Disk " + metric.Name() + " " + std::to_string(stencil.nTh) + "." + std::to_string(stencil.nPh)
              + " " + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z"
              + " s" + std::to_string(sigma) + " Leb" + std::to_string(lebedevStencil.nOrder) + " t" + std::to_string(simTime)
              + ((diskIsHomogeneous) ? " homogeneous" : " inhomogeneous"),
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

// Note:
// -Schwarzschild InsideBH has been modified to 2.4.

// TODO:
// -fix kerr metric initial direction
// -try higher cfl condition
// -check out itp servers
int main()
{
    // Tests:
    // Test_TensorTypes();
    // Test_Grid();
    // Test_Interpolation();
    // Test_Metric();
    // Test_GeodesicEquationSolver();
    // Test_Stencil();
    // Test_LebedevStencil();
    // Test_Quadrature();
    // Test_QuadratureLebedev();
    // Test_QuadratureGaussLegendre();
    // Test_SphericalHarmonicsExpansion();
    // Test_QuaternionRotation();
    // Test_IntensityAt();
    // Test_PhotonFourVelocity();
    // Test_MyAtan2();
    // Test_MySin_MyCos();
    // Test_Camera();
    // Test_Emission( 25, 45, 50,  20,    15,10);

    // Visualization data:
    // BlackHoleCsv();

    // Runs:
    // StreamFlatBeam();
    // StreamFlatStaticSphereWave();
    // StreamFlatDynamicSphereWave();

  //StreamCurvedBeam( nx, ny, nz, nTh, sigma,simTime);
    // StreamCurvedBeam( 26, 46, 51,  15,    15,10);
    // StreamCurvedBeam( 26, 46, 51,  19,   20,10);

    // StreamCurvedBeam( 51, 91, 101,  15,    15,10);
    // StreamCurvedBeam( 51, 91, 101,  19,    20,10);



  //ThinDisk( nx, ny, nz, nTh, sigma,simTime, diskHomogeneous?);
    // Low res, low dir:
    // ThinDisk(106,136, 76,  15,     1,80, true);
    // ThinDisk(106,136, 76,  15,     1,80, false); // done

    // Low res, high dir:
    // ThinDisk(106,136, 76,  19,     1,80, true);
    // ThinDisk(106,136, 76,  19,     1,80, false); // done

    // Middle res, low dir:
    // ThinDisk(148,190,106,  15,     1,80, true);
    // ThinDisk(148,190,106,  15,     1,80, false); // done

    // Middle res, high dir:
    // ThinDisk(148,190,106,  19,     1,80, true);
    ThinDisk(148,190,106,  19,     1,80, false);

    // high res, low dir:
    // ThinDisk(211,271, 151,  15,     1,80, true);
    // ThinDisk(211,271, 151,  15,     1,80, false);

    // high res, high dir:
    // ThinDisk(211,271, 151,  19,     1,80, true);
    // ThinDisk(211,271, 151,  19,     1,80, false);
}