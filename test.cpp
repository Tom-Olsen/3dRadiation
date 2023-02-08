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
    int nx, ny, nz;
    nx = ny = nz = 50;
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
    double data0[grid.nxyz];
    double data1[grid.nxyz];
    double data2[grid.nxyz];
    double data3[grid.nxyz];
    for(int i=0; i<grid.nx; i++)
    for(int j=0; j<grid.ny; j++)
    for(int k=0; k<grid.nz; k++)
    {
        int ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        data0[ijk] = (dist(rOrb,xyz) < 0.2) ? 1 : 0;
        data1[ijk] = (dist(gOrb,xyz) < 0.2) ? 1 : 0;
        data2[ijk] = (dist(bOrb,xyz) < 0.2) ? 1 : 0;
        data3[ijk] = 1;
    }
    grid.WriteFrametoJson(0,data0,data1,data2,data3,0,"output","Test_Grid");
    grid.WriteFrametoCsv (0,data0,data1,data2,data3,0,"output","Test_Grid");
}



void Test_Interpolation()
{
    int nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    double data0[grid.nxyz];
    double data1[grid.nxyz];
    double data2[grid.nxyz];
    double data3[grid.nxyz];
    for(int i=0; i<grid.nx; i++)
    for(int j=0; j<grid.ny; j++)
    for(int k=0; k<grid.nz; k++)
    {
        int ijk = grid.Index(i,j,k);
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
    int nx, ny, nz;
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
    fileOut << "#x, y, z \n";

    int nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-3,-3,-3);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    cout << "Initialization complete." << endl;

    int n = 10;
    for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
    {
        Coord x(start[1] + (i+0.5) * (end[1] - start[1]) / n, start[2] + (j+0.5) * (end[2] - start[2]) / n, start[3]);
        Tensor4x4 g_ll = metric.GetMetric_ll(x);
        Tensor3x3 gamma_ll = metric.GetGamma_ll(x);
        Tensor3x3 delta_ll(1,0,0, 0,1,0, 0,0,1);

        Tensor4 uLF(1,0,0,1);
        uLF = NullNormalize(uLF, g_ll);
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF,x,metric);

        fileOut << x[1] << ", " << x[2] << ", " << x[3] << "\n";
        while(true)
        {
            if (grid.OutsideDomain(x) || metric.InsideBH(x))
                break;
            RK45_GeodesicEquation(10 * grid.dt, x, vLF, metric);
            fileOut << x[1] << ", " << x[2] << ", " << x[3] << "\n";
        }
        cout << "Photon(" << i << "," << j << ") complete." << endl;
    }

    fileOut.close();
}



void Test_Stencil()
{
    Stencil stencil(5,10);
    stencil.Print();

    {
        double d0 = 2.3;
        double d1 = 6.7;
        double x = stencil.Cx(d0,d1);
        double y = stencil.Cy(d0,d1);
        double z = stencil.Cz(d0,d1);
        double theta = stencil.Theta(d0,d1);
        double phi = stencil.Phi(d0,d1);
        cout << theta << ", " << stencil.Theta(x,y,z) << endl;
        cout << phi << ", " << stencil.Phi(x,y,z) << endl;
    }
    cout << endl;
    {
        double theta = 2.1;
        double phi = 4.2;
        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);
        cout << theta << ", " << stencil.Theta(x,y,z) << endl;
        cout << phi << ", " << stencil.Phi(x,y,z) << endl;
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
    int nCoefficients = 1 + 3 + 5;
    int nDir0 = 10;
    int nDir1 = 20;
    Stencil stencil(nDir0,nDir1);

    double data[nDir0 * nDir1];
    for(int d1=0; d1<nDir1; d1++)
    for(int d0=0; d0<nDir0; d0++)
    {
        double x = stencil.Cx(d0,d1);
        double y = stencil.Cy(d0,d1);
        double z = stencil.Cz(d0,d1);

        int d = stencil.Index(d0,d1);
        data[d] = SphericalHarmonics::Y00(x,y,z);

        data[d] += 0.1 * SphericalHarmonics::Y1m1(x,y,z);
        data[d] += 0.2 * SphericalHarmonics::Y10(x,y,z);
        data[d] += 0.3 * SphericalHarmonics::Y1p1(x,y,z);
        
        data[d] += 1.1 * SphericalHarmonics::Y2m2(x,y,z);
        data[d] += 1.2 * SphericalHarmonics::Y2m1(x,y,z);
        data[d] += 1.3 * SphericalHarmonics::Y20(x,y,z);
        data[d] += 1.4 * SphericalHarmonics::Y2p1(x,y,z);
        data[d] += 1.5 * SphericalHarmonics::Y2p2(x,y,z);
    }
    vector<double>c = SphericalHarmonics::GetCoefficients(stencil,data,nCoefficients);
    
    for(int i=0; i<c.size(); i++)
        cout << "c" << i << ": " << c[i] << endl;
    cout << endl;

    cout << "Real Value, Compressed Value:" << endl;
    int d0 = 5;
    int d1 = 15;
    double theta = stencil.Theta(d0,d1);
    double phi = stencil.Phi(d0,d1);
    cout << data[stencil.Index(d0,d1)] << ", " << SphericalHarmonics::GetValue(theta,phi,c,nCoefficients) << endl << endl;
}
void Test_QuadratureLebedev()
{
    LebedevStencil7 lebedev;
    int nDir = lebedev.nDir;

    double data[nDir];
    for(int d=0; d<nDir; d++)
    {
        double x = lebedev.Cx(d);
        double y = lebedev.Cy(d);
        double z = lebedev.Cz(d);

        data[d] = SphericalHarmonics::Y00(x,y,z);

        data[d] += 0.1 * SphericalHarmonics::Y1m1(x,y,z);
        data[d] += 0.2 * SphericalHarmonics::Y10 (x,y,z);
        data[d] += 0.3 * SphericalHarmonics::Y1p1(x,y,z);
        
        data[d] += 1.1 * SphericalHarmonics::Y2m2(x,y,z);
        data[d] += 1.2 * SphericalHarmonics::Y2m1(x,y,z);
        data[d] += 1.3 * SphericalHarmonics::Y20 (x,y,z);
        data[d] += 1.4 * SphericalHarmonics::Y2p1(x,y,z);
        data[d] += 1.5 * SphericalHarmonics::Y2p2(x,y,z);
        
        data[d] += 2.1 * SphericalHarmonics::Y3m3(x,y,z);
        data[d] += 2.2 * SphericalHarmonics::Y3m2(x,y,z);
        data[d] += 2.3 * SphericalHarmonics::Y3m1(x,y,z);
        data[d] += 2.4 * SphericalHarmonics::Y30 (x,y,z);
        data[d] += 2.5 * SphericalHarmonics::Y3p1(x,y,z);
        data[d] += 2.6 * SphericalHarmonics::Y3p2(x,y,z);
        data[d] += 2.7 * SphericalHarmonics::Y3p3(x,y,z);
    }
    vector<double> c = SphericalHarmonics::GetCoefficients(lebedev,data,lebedev.nCoefficients);

    for(int i=0; i<c.size(); i++)
        cout << "c" << i << ": " << c[i] << endl;
    cout << endl;
    
    cout << "Real Value, Compressed Value:" << endl;
    int d = lebedev.nDir / 2;
    double theta = lebedev.Theta(d);
    double phi = lebedev.Phi(d);
    cout << data[d] << ", " << SphericalHarmonics::GetValue(theta,phi,c,lebedev.nCoefficients) << endl << endl;

    cout << "         nDir = " << lebedev.nDir << endl;
    cout << "       nOrder = " << lebedev.nOrder << endl;
    cout << "nCoefficients = " << lebedev.nCoefficients << endl;
}



void Test_SphericalHarmonicsExpansion()
{
    int nx, ny, nz;
    nx = ny = nz = 20;
    Coord start(0,0,0);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    Stencil stencil(15,20);
    LebedevStencil3 lebedevStencil;

    Radiation radiation(metric, stencil, lebedevStencil, StreamingType::CurvedDynamic);
    radiation.UpdateSphericalHarmonicsCoefficients();
    
    ofstream file("output/Test_SphericalHarmonicsExpansion.csv");
    file << "#x, y, z, s \n";
	for(int k=0; k<grid.nz; k+=2)
	for(int j=0; j<grid.ny; j+=2)
	for(int i=0; i<grid.nx; i+=2)
    {
		if(i == grid.nx/2 && j == grid.ny/2 && k == grid.nz/2)
        {
	        std::vector<double> cS(lebedevStencil.nCoefficients);
	        for(int f=0; f<lebedevStencil.nCoefficients; f++)
	        	cS[f] = radiation.coefficientsS[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
                
    		for(int i=0; i<cS.size(); i++)
    		    std::cout << "cS" << i << ": " << cS[i] << std::endl;
    		std::cout << std::endl;
        }

        for(int d1=0; d1<stencil.nPh; d1++)
        for(int d0=0; d0<stencil.nTh; d0++)
        {
            double theta = stencil.Theta(d0,d1);
            double phi = stencil.Phi(d0,d1);
            double s = radiation.GetFrequencyShift(i,j,k,theta,phi);
            Coord xyz = radiation.GetTempCoordinate(i,j,k,theta,phi);
            Tensor3 v = radiation.GetTemp3Velocity(i,j,k,theta,phi);
            if(!metric.InsideBH(xyz))
                file << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << s << "\n";
        }
        // North Pole:
        {
            double theta = 0;
            double phi = 0;
            double s = radiation.GetFrequencyShift(i,j,k,theta,phi);
            Coord xyz = radiation.GetTempCoordinate(i,j,k,theta,phi);
            Tensor3 v = radiation.GetTemp3Velocity(i,j,k,theta,phi);
            if(!metric.InsideBH(xyz))
                file << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << s << "\n";
        }
        // South Pole:
        {
            double theta = M_PI;
            double phi = 0;
            double s = radiation.GetFrequencyShift(i,j,k,theta,phi);
            Coord xyz = radiation.GetTempCoordinate(i,j,k,theta,phi);
            Tensor3 v = radiation.GetTemp3Velocity(i,j,k,theta,phi);
            if(!metric.InsideBH(xyz))
                file << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << s << "\n";
        }
    }
    file.close();
}



void Test_StreamFlatStatic()
{
    // Create Radiation object:
    int nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil stencil(6,10);
    LebedevStencil3 lebedevStencil;
    Radiation radiation(metric, stencil, lebedevStencil, StreamingType::FlatStatic);

    // Initial Data:
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double r = xyz.Radius();
        if (r < 0.2)
        {
            radiation.isInitialGridPoint[ijk] = true;
            radiation.initialE[ijk] = 1;
        }
    }

    // Start simulation:
    Config config =
    {
        .name = "Test_StreamFlatStatic",
        .simTime = 1,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true
    };
    radiation.RunSimulation(config);
}



void Test_QuaternionRotation()
{
    using namespace glm;
    vec3 from(0,0,1);
    vec3 to(1,0,0);
    quat q = quat(from, to);
    // quat qInv = glm::conjugate(q);
    quat qInv = quat(to, from);

    Stencil stencil(30,40);

    ofstream file0("output/Test_QuaternionRotation Normal");
    ofstream file1("output/Test_QuaternionRotation Rotated");
    file0 << "#x, y, z" << endl;
    file1 << "#x, y, z" << endl;
    
    for(int d1=0; d1<stencil.nPh; d1++)
    for(int d0=0; d0<stencil.nTh; d0++)
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
    // This test is simulating the Radiation::IntensityAt(int ijk, Tensor3 vTempIF) function.
    Stencil stencil(5,8);
    LebedevStencil11 lebedev;

    double I[stencil.nThPh];
    double Inorth;
    double Isouth;

    ofstream file0("output/Test_IntensityAtStencil");
    ofstream file1("output/Test_IntensityAtLebedev");
    file0 << "x, y, z, color\n";
    file1 << "x, y, z, color\n";

    for(int d1=0; d1<stencil.nPh; d1++)
    for(int d0=0; d0<stencil.nTh; d0++)
    {
        int d = stencil.Index(d0,d1);
        Tensor3 xyz = stencil.Cxyz(d0,d1);
        I[d] = 1.0 + 0.75 * SphericalHarmonics::Y1p1(xyz) + 0.75 * SphericalHarmonics::Y10 (xyz) + 0.75 * SphericalHarmonics::Y1m1(xyz);
        file0 << I[d] * xyz[1] << ", " << I[d] * xyz[2] << ", " << I[d] * xyz[3] << ", " << 1 << "\n";
    }
    Inorth = 1.0 + 0.75 * SphericalHarmonics::Y1p1(0, 0, 1) + 0.75 * SphericalHarmonics::Y10 (0, 0, 1) + 0.75 * SphericalHarmonics::Y1m1(0, 0, 1);
    Isouth = 1.0 + 0.75 * SphericalHarmonics::Y1p1(0, 0,-1) + 0.75 * SphericalHarmonics::Y10 (0, 0,-1) + 0.75 * SphericalHarmonics::Y1m1(0, 0,-1);
    file0 << Inorth * 0 << ", " << Inorth * 0 << ", " << Inorth * 1 << ", "  << 0 << "\n";
    file0 << Isouth * 0 << ", " << Isouth * 0 << ", " << Isouth * -1 << ", " << 0 << "\n";

    for(int d=0; d<lebedev.nDir; d++)
    {
        Tensor3 xyz = lebedev.Cxyz(d);
	    double t = stencil.d0(xyz.Theta());
	    double p = stencil.d1(xyz.Phi());
        
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
void Test_StreamFlatDynamic()
{
    // Create Radiation object:
    int nx, ny, nz;
    nx = ny = nz = 100;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil stencil(6,10);
    LebedevStencil3 lebedevStencil;
    Radiation radiation(metric, stencil, lebedevStencil, StreamingType::FlatDynamic);

    // Initial Data:
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
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
        .name = "Test_StreamFlatDynamic",
        .simTime = 1,
        .writeFrequency = 5,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true
    };
    radiation.RunSimulation(config);
}



int main()
{
    // Test_TensorTypes();
    // Test_Grid();
    // Test_Interpolation();
    // Test_Metric();
    // Test_GeodesicEquationSolver();
    // Test_Stencil();
    // Test_LebedevStencil();
    // Test_Quadrature();
    // Test_QuadratureLebedev();
    // Test_SphericalHarmonicsExpansion();
    // Test_StreamFlatStatic();
    // Test_QuaternionRotation();
    // Test_IntensityAt();
    Test_StreamFlatDynamic();
}