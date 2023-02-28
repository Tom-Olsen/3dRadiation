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
    // grid.WriteFrametoJson(0,data0,data1,data2,data3,0,"output","Test_Grid");
    grid.WriteFrametoCsv (0,data0,data1,data2,data3,0,"output","Test_Grid");
}



void Test_Interpolation()
{
    int nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    RealBuffer data0(grid.nxyz);
    RealBuffer data1(grid.nxyz);
    RealBuffer data2(grid.nxyz);
    RealBuffer data3(grid.nxyz);
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
    fileOut << "#x, y, z, s \n";

    int nx, ny, nz;
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
        Tensor3 cxyz = stencil.Cxyz(d0,d1);
        cout << theta << ", " << stencil.Theta(x,y,z) << ", " << cxyz.Theta() << endl;
        cout << phi << ", " << stencil.Phi(x,y,z) << ", " << cxyz.Phi() << endl;
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
	for(int k=2; k<grid.nz-2; k+=2)
	for(int j=2; j<grid.ny-2; j+=2)
	for(int i=2; i<grid.nx-2; i+=2)
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

        int ijk = grid.Index(i,j,k);
        double alpha = metric.GetAlpha(ijk);
        Coord center = grid.xyz(i,j,k);
        for(int d1=0; d1<stencil.nPh; d1++)
        for(int d0=0; d0<stencil.nTh; d0++)
        {
            // Get pos and vel by Spherical Harmonic Expansion:
            {
                double theta = stencil.Theta(d0,d1);
                double phi = stencil.Phi(d0,d1);
                double s = radiation.GetFrequencyShift(i,j,k,theta,phi);
                Coord xyz = radiation.GetTempCoordinate(i,j,k,theta,phi);
                Tensor3 v = radiation.GetTemp3Velocity(i,j,k,theta,phi);
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



void Test_StreamFlatStaticSphereWave()
{
    // Create Radiation object:
    int nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil stencil(6,10);
    stencil.sigma = 1.0;    // not used in this test due to direction: nx=ny=nz=0
    LebedevStencil3 lebedevStencil;
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);

    // Initial Data:
    PARALLEL_FOR(3)
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
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
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



void Test_StreamFlatDynamicSphereWave()
{
    // Create Radiation object:
    int nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-1,-1,-1);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil stencil(6,10);
    stencil.sigma = 1.0;    // not used in this test due to direction: nx=ny=nz=0
    LebedevStencil3 lebedevStencil;
    Camera camera;
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);

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
            radiation.initialNx[ijk] = 0;
            radiation.initialNy[ijk] = 0;
            radiation.initialNz[ijk] = 0;
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
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



void Test_StreamFlatBeam()
{
    // Create Radiation object:
    int nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(-2,0,0);
    Coord end(2,4,4);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil stencil(15,30);
    stencil.sigma = 100;
    LebedevStencil3 lebedevStencil;
    Camera camera;
    // Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatDynamic);
    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::FlatStatic);

    // Initial Data:
    double beamWidth = 0.5;
    Coord beamCenter(0,3,0);
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
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
            radiation.initialKappa0[ijk] = 0;
            radiation.initialKappa1[ijk] = 0;
            radiation.initialKappaA[ijk] = 0;
            radiation.initialEta[ijk] = 0;
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



void Test_StreamCurvedBeam(int nx, int ny, int nz, int nTh, int nPh, int sigma, int simTime)
{
    // Grid, Metric, Stencil:
    Coord start(-1,0,-0.1);
    Coord end(1,3.6,3.9);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);
    Stencil stencil(nTh,nPh);
    stencil.sigma = sigma;
    LebedevStencil5 lebedevStencil;

    // Camera:
    int resX = 100;
    int resY = 200;
    int width = 2;
    int height = 4;
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
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        double x = xyz[1];
        double y = xyz[2];
        double z = xyz[3];
        if (-0.25 < x && x < 0.25
          && 3.00 < y && y < 3.50
          && z < 0.00)
        {
            double alpha = metric.GetAlpha(ijk);
            Tensor4 u(alpha,0,0,alpha);
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



void Test_Camera()
{
    // Black Hole with Thin Disk Camera Settings:
    {
        int resX = 400;
        int resY = 200;
        int width = 28;
        int height = 14;

        Coord position(0,20,3);
        double degreeToRadians = 2.0 * M_PI / 360.0;
        double angleX = 100 * degreeToRadians;
        double angleY = 0 * degreeToRadians;
        double angleZ = 0 * degreeToRadians;
        glm::vec3 eulerAngles(angleX,angleY,angleZ);

        Camera camera(resX, resY, width, height, position, eulerAngles);

        for(int ij=0; ij<camera.pixelCount; ij++)
        {
            int i = ij % camera.resX;
            int j = ij / camera.resX;
            camera.image[ij] = i / (resX-1.0) * j / (resY-1.0);
        }

        camera.WriteImagetoCsv(0, 0, "output");
    }
    
    // Curved Beam Close Settings:
    {
        int resX = 100;
        int resY = 50;
        int width = 2;
        int height = 1;

        Coord position(0,3/sqrt(2.0),3/sqrt(2.0));
        double degreeToRadians = 2.0 * M_PI / 360.0;
        double angleX = -145 * degreeToRadians;
        double angleY = 0 * degreeToRadians;
        double angleZ = 0 * degreeToRadians;
        glm::vec3 eulerAngles(angleX,angleY,angleZ);

        Camera camera(resX, resY, width, height, position, eulerAngles);

        for(int ij=0; ij<camera.pixelCount; ij++)
        {
            int i = ij % camera.resX;
            int j = ij / camera.resX;
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

    ofstream file0("output/BlackHole.csv");
    ofstream file1("output/ThinDisk.csv");
    file0 << "#x,y,z,value\n";
    file1 << "#x,y,z,value\n";

    // Add one point with different value so min != max.
    file0 << "0,0,0,1" << "\n";
    file1 << "0,0,0,0" << "\n";

    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
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



void Test_Emission(int nx, int ny, int nz, int nTh, int nPh, int sigma, int simTime)
{
    // Grid, Metric, Stencil:
    Coord start(-1,0,-0.1);
    Coord end(1,3.6,3.9);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);
    Stencil stencil(nTh,nPh);
    stencil.sigma = sigma;
    LebedevStencil5 lebedevStencil;

    // Camera:
    int resX = 100;
    int resY = 200;
    int width = 2;
    int height = 4;
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
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
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



void Test_ThinDisk(int nx, int ny, int nz, int nTh, int nPh, int sigma, int simTime)
{
    // Black Hole and Thin Disk:
    double m = 1;
    double a = 0;
    double r = 2 * m;
    double diskInner = 3 * r;   //  6
    double diskOuter = 6 * r;   // 12

    // Grid, Metric, Stencil:
    Coord start(-14,-14,-10);
    Coord   end( 14, 22, 10);
    Grid grid(nx, ny, nz, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, m, a);
    Stencil stencil(nTh,nPh);
    stencil.sigma = sigma;
    LebedevStencil5 lebedevStencil;

    // Camera:
    int resX = 400;
    int resY = 200;
    int width = 26;
    int height = 13;
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
    for(int k=0; k<grid.nz; k++)
    for(int j=0; j<grid.ny; j++)
    for(int i=0; i<grid.nx; i++)
    {
        int ijk = grid.Index(i,j,k);
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
            radiation.initialEta[ijk] = 1;
        }
    }

    // Start simulation:
    Config config =
    {
        .name = "Thin Disk " + metric.Name() + " " + std::to_string(stencil.nTh) + "." + std::to_string(stencil.nPh)
              + " s" + std::to_string(sigma) + " Leb" + std::to_string(lebedevStencil.nOrder) + " t" + std::to_string(simTime),
        .simTime = (double)simTime,
        .writeFrequency = 20,
        .updateSphericalHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
        .useCamera = true
    };
    cout << "starting run" << endl;
    radiation.RunSimulation(config);
}


// Note:
// -Schwarzschild InsideBH has been modified to 2.4.
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
    // Test_StreamFlatStaticSphereWave();
    // Test_QuaternionRotation();
    // Test_IntensityAt();
    // Test_StreamFlatDynamicSphereWave();
    // Test_StreamFlatBeam();

  //Test_StreamCurvedBeam( nx, ny, nz, nTh,nPh, sigma,simTime);
    Test_StreamCurvedBeam( 25, 45, 50,  20, 40,    15,10);
    // Test_StreamCurvedBeam( 50, 90,100,  20, 40,    15,10);

    // Test_Emission( 25, 45, 50,  20, 40,    15,10);

    // Test_Camera();
    // BlackHoleCsv();
  //Test_ThinDisk( nx, ny, nz, nTh,nPh, sigma,simTime);
    // Test_ThinDisk(106,136, 76,  20, 40,     1,200);
}