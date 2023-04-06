#include <iostream>
#include "src/Includes.hh"

using namespace std;






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



void Test_PhotonVelocity()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(1,1,1);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    KerrSchild metric(grid, 1.0, 0.0);
    
    Coord xyz(1.5, 2, 2.5);
    double alpha = metric.GetAlpha(xyz);
    Tensor4x4 g_ll = metric.GetMetric_ll(xyz);
    Tensor3x3 eta3x3_ll = metric.GetMinkowskiGamma_ll(xyz);
    Tensor4x4 eta4x4_ll = metric.GetMinkowskiMetric_ll(xyz);
    Tensor3x3 gamma_ll = metric.GetGamma_ll(xyz);

    // Four velocity defined in LF and transformed to three velocity in IF:
    {
        Tensor4 uLF(1,0.2,0.4,0.6);
        uLF = NullNormalize(uLF,g_ll);
        Tensor3 vIF = Vec3ObservedByEulObs<LF,IF>(uLF,xyz,metric);

        uLF.Print(" uLF ");
        PrintDouble(Norm2(uLF,g_ll),"|uLF|");
        vIF.Print(" vIF ");
        PrintDouble(Norm2(vIF,eta3x3_ll),"|vIF|");
    }cout << endl;
    // Four velocity defined in IF and transformed to three velocity in LF:
    {
        Tensor4 uIF(1*alpha,0.2*alpha,0.4*alpha,0.6*alpha);
        Tensor3 vLF = Vec3ObservedByEulObs<IF,LF>(uIF, xyz, metric);

        uIF.Print(" uIF ");
        PrintDouble(Norm2(uIF,eta4x4_ll),"|uIF|");
        vLF.Print(" vLF ");
        PrintDouble(Norm2(vLF,gamma_ll),"|vLF|");
    }cout << endl;
    // Four velocity defined in LF and transformed to three velocity in LF:
    {
        Tensor4 uLF(1,0.2,0.4,0.6);
        uLF = NullNormalize(uLF,g_ll);
        Tensor3 vLF = Vec3ObservedByEulObs<LF,LF>(uLF, xyz, metric);

        uLF.Print(" uLF ");
        PrintDouble(Norm2(uLF,g_ll),"|uLF|");
        vLF.Print(" vLF ");
        PrintDouble(Norm2(vLF,gamma_ll),"|vLF|");
    }
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



void Test_Stencil(int nOrder)
{
    MyStencil stencil(nOrder);

    double d0 = 2.3;
    double d1 = 6.7;
    double x = stencil.Cx(d0,d1);
    double y = stencil.Cy(d0,d1);
    double z = stencil.Cz(d0,d1);
    double theta = stencil.Theta(d0,d1);
    double phi = stencil.Phi(d0,d1);
    Tensor3 cxyz = stencil.C(d0,d1);
    cout << theta << ", " << cxyz.Theta() << endl;
    cout << phi << ", " << cxyz.Phi() << endl;
    cout << endl;

    LebedevStencil lebStencil(nOrder);
    for(int d=0; d<lebStencil.nDir; d++)
    {
        // cout << Format(lebStencil.Theta(d)) << ", " << Format(lebStencil.Phi(d)) << ", " << Format(lebStencil.W(d)) << ", " << endl;
        // cout << Format(lebStencil.Cx(d)) << ", " << Format(lebStencil.Cy(d)) << ", " << Format(lebStencil.Cz(d)) << ", " << Format(lebStencil.W(d)) << endl;
        // cout << Format(lebStencil.Cx(d)) << ", " << Format(lebStencil.Cy(d)) << ", " << Format(lebStencil.Cz(d)) << endl;
    }
    cout << endl;
    cout << endl;

    //ofstream file("output/sphere.txt");
    //file << "x,y,z,c\n";
    //int n = 500;
    //double minMaxDot = 1;
    //for(int i=0; i<n; i++)
    //for(int j=0; j<2*n; j++)
    //{
    //    double theta = M_PI * (i + 0.5) / (double)n;
    //    double phi = 2.0 * M_PI * j / (double)n;
    //    Tensor3 c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    //    double maxDot = -1;
    //    for(int d=0; d<lebStencil.nDir; d++)
    //    {
    //        double dot = Tensor3::Dot(lebStencil.C(d), c);
    //        maxDot = std::max(maxDot,dot);
    //    }
    //    file << c[1] << "," << c[2] << "," << c[3] << "," << maxDot << "\n";
    //    minMaxDot = std::min(minMaxDot,maxDot);
    //}
    //file.close();

    lebStencil.connectedTriangles.Print();
    cout << endl;
    lebStencil.connectedVertices.Print();
    cout << endl;
    cout << "minMaxDot: " << lebStencil.minMaxDot << endl;
}



void Test_TriangleIntersection()
{
    Tensor3 v0( 0, 1, 1);
    Tensor3 v1(-1,-1, 1);
    Tensor3 v2( 1,-1, 1);

    Tensor3 color0(1,0,0);
    Tensor3 color1(0,1,0);
    Tensor3 color2(0,0,1);

    Tensor3 origin(0,0,0);
    Tensor3 intersection;
    Tensor3 weights;
    int n = 101;

    ofstream file("output/barycentricWeights.csv");
    file << "x,y,z,r,g,b\n";
    for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    {
        double x = 2.0 * i / (n - 1.0) - 1.0;
        double y = 2.0 * j / (n - 1.0) - 1.0;
        Tensor3 direction(x,y,1);
        //double norm = direction.EuklNorm();
        //direction[1] /= norm;
        //direction[2] /= norm;
        //direction[3] /= norm;
        // if(RayTriangleIntersection(origin, direction, v0, v1, v2, intersection))
            // cout << intersection[1] << "," << intersection[2] << "," << intersection[3] << endl;
        RayTriangleIntersection(origin, direction, v0, v1, v2, intersection);
        if(BarycentricWeights(origin, direction, v0, v1, v2, weights))
        {
            Tensor3 color = color0 * weights[1] + color1 * weights[2] + color2 * weights[3];
            file << intersection[1] << "," << intersection[2] << "," << intersection[3] << ","
                 << color[1] << "," << color[2] << "," << color[3] << "\n";
        }
    }
    file.close();
}



void Test_SphericalBarycentricWeights()
{
    Tensor3 a(1,0,0);
    Tensor3 b(0,1,0);
    Tensor3 c(0,0,1);
    
    Tensor3 p(1,1,1);
    double norm = p.EuklNorm();
    p[1] /= norm;
    p[2] /= norm;
    p[3] /= norm;

    Tensor3 weights;
    SphericalBarycentricWeights(p,a,b,c,weights);
}



void Test_Quadrature(const Stencil& stencil)
{
    cout << "nDir = " << stencil.nDir << endl;
    cout << "nOrder = " << stencil.nOrder << endl;
    cout << "nCoefficients = " << stencil.nCoefficients << endl;

    double data[stencil.nDir];
    for(size_t d=0; d<stencil.nDir; d++)
    {
        Tensor3 dir = stencil.C(d);
        data[d] = 0;
        for(int k=0; k<stencil.nCoefficients; k++)
            data[d] += (k+1) * SphericalHarmonicsXyz::Y(k,dir);
    }

    double coefficients[stencil.nDir];
    SphericalHarmonicsXyz::GetCoefficients(stencil,data,coefficients);
    // vector<double> coefficients = SphericalHarmonicsXyz::GetCoefficients(*stencil,data);
    
    for(size_t i=0; i<stencil.nCoefficients; i++)
        cout << "c" << i << ": " << coefficients[i] << endl;
    cout << endl;

    cout << "Real Value, Compressed Value:" << endl;
    size_t d = stencil.nDir/2;
    cout << data[d] << ", " << SphericalHarmonicsXyz::GetValue(stencil.C(d),coefficients,stencil.nCoefficients) << endl << endl;
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
    MyStencil stencil(15);
    LebedevStencil lebedevStencil(5);
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
		//if(i == grid.nx/2 && j == grid.ny/2 && k == grid.nz/2)
        //{
	    //    std::vector<double> cS(lebedevStencil.nCoefficients);
	    //    for(size_t f=0; f<lebedevStencil.nCoefficients; f++)
	    //    	cS[f] = radiation.coefficientsS[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
        //        
    	//	for(size_t i=0; i<cS.size(); i++)
    	//	    std::cout << "cS" << i << ": " << cS[i] << std::endl;
    	//	std::cout << std::endl;
        //}

        size_t ijk = grid.Index(i,j,k);
        double alpha = metric.GetAlpha(ijk);
        Coord center = grid.xyz(i,j,k);
        for(size_t d1=0; d1<stencil.nPh; d1++)
        for(size_t d0=0; d0<stencil.nTh; d0++)
        {
            // Get pos and vel by Spherical Harmonic Expansion:
            {
                Tensor3 direction = stencil.C(d0,d1);
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
                Tensor3 c = stencil.C(d0,d1);
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
}



void Test_QuaternionRotation()
{
    using namespace glm;
    vec3 from(0,0,1);
    vec3 to(1,0,0);
    quat q = quat(from, to);
    quat qInv = quat(to, from);

    MyStencil stencil(15);

    ofstream file0("output/Test_QuaternionRotation Normal.csv");
    ofstream file1("output/Test_QuaternionRotation Rotated.csv");
    file0 << "#x, y, z" << endl;
    file1 << "#x, y, z" << endl;
    
    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        Tensor3 p = stencil.C(d0,d1);
        file0 << p[1] << ", " << p[2] << ", " << p[3] << endl;
        p = q * p;
        file1 << p[1] << ", " << p[2] << ", " << p[3] << endl;
    }

    file0.close();
    file1.close();

    Tensor3 p(1,2,3);
    Tensor3 pRot = q * p;
    Tensor3 pRotInvRot = qInv * pRot;
    Tensor3 pRotInvRot2 = Invert(q) * pRot;

    p.Print("   p");
    pRot.Print("pRot");
    pRotInvRot.Print("  p?");
    pRotInvRot2.Print("  p?");
}



void Test_IntensityAt()
{
    // This test is simulating the Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF) function.
    MyStencil stencil(5);
    LebedevStencil lebedev(11);

    double I[stencil.nDir];
    double Inorth;
    double Isouth;

    ofstream file0("output/Test_IntensityAtStencil.csv");
    ofstream file1("output/Test_IntensityAtLebedev.csv");
    file0 << "#x, y, z, color\n";
    file1 << "#x, y, z, color\n";

    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {
        size_t d = stencil.Index(d0,d1);
        Tensor3 dir = stencil.C(d0,d1);
        I[d] = 1.0 + 0.75 * SphericalHarmonicsXyz::Y(1,dir) + 0.75 * SphericalHarmonicsXyz::Y(2,dir) + 0.75 * SphericalHarmonicsXyz::Y(3,dir);
        file0 << I[d] * dir[1] << ", " << I[d] * dir[2] << ", " << I[d] * dir[3] << ", " << 1 << "\n";
    }
    Tensor3 dir(0,0,1);
    Inorth = 1.0 + 0.75 * SphericalHarmonicsXyz::Y(1,dir) + 0.75 * SphericalHarmonicsXyz::Y(2,dir) + 0.75 * SphericalHarmonicsXyz::Y(3,dir);
    dir = Tensor3(0,0,-1);
    Isouth = 1.0 + 0.75 * SphericalHarmonicsXyz::Y(1,dir) + 0.75 * SphericalHarmonicsXyz::Y(2,dir) + 0.75 * SphericalHarmonicsXyz::Y(3,dir);
    file0 << Inorth * 0 << ", " << Inorth * 0 << ", " << Inorth * 1 << ", "  << 0 << "\n";
    file0 << Isouth * 0 << ", " << Isouth * 0 << ", " << Isouth * -1 << ", " << 0 << "\n";

    for(size_t d=0; d<lebedev.nDir; d++)
    {
        Tensor3 dir = lebedev.C(d);
	    double t = stencil.d0(dir.Theta());
	    double p = stencil.d1(dir.Phi());
        
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
	    double d00 = stencil.d0(0);
        if (t0 == -1)			// North
		    tfrac = 1.0 / abs(d00) * (tfrac - (1.0 + d00));
        if (t1 == stencil.nTh)	// South
		tfrac /= abs(d00);

	    double I = BilinearInterpolation(tfrac, pfrac, I00, I01, I10, I11);
        double Iexpected = 1.0 + 0.75 * SphericalHarmonicsXyz::Y(1,dir) + 0.75 * SphericalHarmonicsXyz::Y(2,dir) + 0.75 * SphericalHarmonicsXyz::Y(3,dir);
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
        file1 << I * dir[1] << ", " << I * dir[2] << ", " << I * dir[3] << ", " << color << "\n";
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
		default: { exit_on_error("Invalid IntensityProfile"); }
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



void Test_HarmonicsBenchmark()
{
    LebedevStencil stencil(31);
    double* data1 = new double[stencil.nDir];
    double* data2 = new double[stencil.nDir];
    double* coefficients1 = new double[stencil.nCoefficients];
    double* coefficients2 = new double[stencil.nCoefficients];
    double errorTheshhold = 1e-5;

	Profiler::Session& session = Profiler::Session::Get();
    {
        PROFILE_SCOPE("ThPH");
        {
            PROFILE_SCOPE("Set Data ThPh");
            for(int d=0; d<stencil.nDir; d++)
            {
                float theta = stencil.Theta(d);
                float phi = stencil.Phi(d);

                data1[d] = 0;
                for(int k=0; k<stencil.nCoefficients; k++)
                    data1[d] += (k+1) * SphericalHarmonicsThPh::Y(k,theta,phi);
            }
        }
        {
            PROFILE_SCOPE("Get Coefficients ThPh");
            SphericalHarmonicsThPh::GetCoefficients(stencil, data1, coefficients1);
            for(int k=0; k<stencil.nCoefficients; k++)
                cout << "c[" << k << "]: " << coefficients1[k] << endl;
        }
        {
            PROFILE_SCOPE("Get Values ThPh");
            for(int d=0; d<stencil.nDir; d++)
            {
                float theta = stencil.Theta(d);
                float phi = stencil.Phi(d);
                data2[d] = SphericalHarmonicsThPh::GetValue(theta,phi,coefficients1,stencil.nCoefficients);
            }
        }
        {
            PROFILE_SCOPE("Check Results ThPh");
            int errors = 0;
            for(int d=0; d<stencil.nDir; d++)
            {
                if(abs(data2[d] - data1[d]) > errorTheshhold)
                {
                    // cout << "data1[" << d << "] = " << Format(data1[d]) << ", data2[" << d << "] = " << Format(data2[d]) << ", diff = " << Format(abs(data2[d] - data1[d])) << endl;
                    errors++;
                }
            }
            cout << "errors ThPh: " << errors << endl;
        }
    }
    {
        PROFILE_SCOPE("Xyz");
        {
            PROFILE_SCOPE("Set Data Xyz");
            for(int d=0; d<stencil.nDir; d++)
            {
                float x = stencil.Cx(d);
                float y = stencil.Cy(d);
                float z = stencil.Cz(d);

                data1[d] = 0;
                for(int k=0; k<stencil.nCoefficients; k++)
                    data1[d] += (k+1) * SphericalHarmonicsXyz::Y(k,x,y,z);
            }
        }
        {
            PROFILE_SCOPE("Get Coefficients Xyz");
            SphericalHarmonicsXyz::GetCoefficients(stencil, data1, coefficients2);
            for(int k=0; k<stencil.nCoefficients; k++)
                cout << "c[" << k << "]: " << coefficients2[k] << endl;
        }
        {
            PROFILE_SCOPE("Get Values Xyz");
            for(int d=0; d<stencil.nDir; d++)
            {
                float x = stencil.Cx(d);
                float y = stencil.Cy(d);
                float z = stencil.Cz(d);
                data2[d] = SphericalHarmonicsXyz::GetValue(x,y,z,coefficients2,stencil.nCoefficients);
            }
        }
        {
            PROFILE_SCOPE("Check Results Xyz");
            int errors = 0;
            for(int d=0; d<stencil.nDir; d++)
            {
                if(abs(data2[d] - data1[d]) > errorTheshhold)
                {
                    // cout << "data1[" << d << "] = " << Format(data1[d]) << ", data2[" << d << "] = " << Format(data2[d]) << ", diff = " << Format(abs(data2[d] - data1[d])) << endl;
                    errors++;
                }
            }
            cout << "errors Xyz: " << errors << endl;
        }
    }
	session.End();
    cout << endl;

	std::vector<std::string> names = session.GetAllFunctionNames();
	for(int i=0; i<names.size(); i++)
        session.PrintFunctionDuration(names[i]);
}



void Test_UnstructuredMatrix()
{
    UnstructuredMatrix<double> unstMat;
    double a[3] = {0,1,2};
    double b[4] = {10,11,12,13};
    double c[2] = {20,21};
    double d[5] = {30,31,32,33,34};

    unstMat.AddRow(a,3);
    unstMat.AddRow(b,4);
    unstMat.AddRow(c,2);
    unstMat.AddRow(d,5);

    unstMat.Print();
    unstMat.PrintFlat();

    for(int d=0; d<4; d++)
    {
        size_t start = unstMat.Start(d);
        size_t end = unstMat.End(d);
        cout << d << ": ";
        for(int p=start; p<end; p++)
            cout << unstMat[p] << ",";
        cout << endl;
    }
}



void Test_ConvexHull()
{
    //// Cube test:
    //{
    //    std::vector<Vector3> vertices;
    //    vertices.push_back(Vector3(0,0,0));
    //    vertices.push_back(Vector3(0,0,1));
    //    vertices.push_back(Vector3(0,1,0));
    //    vertices.push_back(Vector3(0,1,1));
    //    vertices.push_back(Vector3(1,0,0));
    //    vertices.push_back(Vector3(1,0,1));
    //    vertices.push_back(Vector3(1,1,0));
    //    vertices.push_back(Vector3(1,1,1));
    //
    //    ConvexHull convexHull(vertices);
    //    Mesh mesh = convexHull.GetMesh();
    //    mesh.WriteToUnityCsv("","cube");
    //}

    // Stencil test:
    // int i = 5;
    for(int i=3; i<=31; i+=2)
    {
        LebedevStencil stencil(i);
        std::vector<Vector3> vertices;
        vertices.reserve(stencil.nDir);

        for(int d=0; d<stencil.nDir; d++)
        {
            Tensor3 c = stencil.C(d);
            Vector3 v(c[1],c[2],c[3]);
            vertices.push_back(v);
        }
        
        ConvexHull convexHull(vertices);
        Mesh mesh = convexHull.GetMesh();
        // mesh.WriteToUnityCsv("output","Lebedev" + to_string(i));

        convexHull.OriginalOrdering(vertices);
        mesh = convexHull.GetMesh();
        mesh.WriteToUnityCsv("output","Lebedev" + to_string(i) + "ordered");
    }
}



void Test_GetTriangle()
{
    LebedevStencil stencil(35);

    // Random direction:
    int N = 100;
    for(int i=0; i<N; i++)
    for(int j=0; j<2*N; j++)
    {
        double theta = M_PI * (i + 0.5) / (double)N;
        double phi = 2.0 * M_PI * j / (double)N;
        Tensor3 n(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

        Vector3Int triangleA;
        Vector3Int triangleB;

        // Find hit triangle with brute force:
        {
            int d;
            double dotMax = -1;
            for(int p=0; p<stencil.nDir; p++)
            {
                double dot = Tensor3::Dot(n,stencil.C(p));
                if(dot > dotMax)
                {
                    dotMax = dot;
                    d = p;
                }
            }

            // Get all triangles connected to direction d:
            int indexStart = stencil.connectedTriangles.Start(d);
            int indexEnd = stencil.connectedTriangles.End(d);


            // Find triangle that gets hit by vTempIF and corresponding barycentric weights:
            Tensor3 weights;
            Tensor3 rayOrigin(0,0,0);

            for(int p=indexStart; p<indexEnd; p++)
            {
                triangleA = stencil.connectedTriangles[p];
                Tensor3 v0 = stencil.C(triangleA[0]);
                Tensor3 v1 = stencil.C(triangleA[1]);
                Tensor3 v2 = stencil.C(triangleA[2]);

                if (BarycentricWeights(rayOrigin, n, v0, v1, v2, weights))
                    break;
            }
        }

        // Find hit triangle with minMaxDot method:
        {
            int d;
            for(int p=0; p<stencil.nDir; p++)
            {
                double dot = Tensor3::Dot(n,stencil.C(p));
                if(dot > stencil.minMaxDot)
                {
                    d = p;
                    break;
                }
            }
            
            int indexStart = stencil.connectedTriangles.Start(d);
            int indexEnd = stencil.connectedTriangles.End(d);
            double dotMax = -1;
            for(int p=indexStart; p<indexEnd; p++)
            {
                Vector3Int triangle = stencil.connectedTriangles[p];
                double dot0 = Tensor3::Dot(n,stencil.C(triangle[0]));
                double dot1 = Tensor3::Dot(n,stencil.C(triangle[1]));
                double dot2 = Tensor3::Dot(n,stencil.C(triangle[2]));
                if(dot0 > dotMax)
                {
                    dotMax = dot0;
                    d = triangle[0];
                }
                if(dot1 > dotMax)
                {
                    dotMax = dot1;
                    d = triangle[1];
                }
                if(dot2 > dotMax)
                {
                    dotMax = dot2;
                    d = triangle[2];
                }
            }

            // Get all triangles connected to direction d:
            indexStart = stencil.connectedTriangles.Start(d);
            indexEnd = stencil.connectedTriangles.End(d);


            // Find triangle that gets hit by vTempIF and corresponding barycentric weights:
            Tensor3 weights;
            Tensor3 rayOrigin(0,0,0);

            for(int p=indexStart; p<indexEnd; p++)
            {
                triangleB = stencil.connectedTriangles[p];
                Tensor3 v0 = stencil.C(triangleB[0]);
                Tensor3 v1 = stencil.C(triangleB[1]);
                Tensor3 v2 = stencil.C(triangleB[2]);

                if (BarycentricWeights(rayOrigin, n, v0, v1, v2, weights))
                    break;
            }
        }
        if(triangleA != triangleB)
        {
            cout << "triangle: " << triangleA << endl;
            cout << "triangle: " << triangleB << endl;
            cout << endl;
        }
    }
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
    // Test_Metric();
    // Test_PhotonVelocity();
    // Test_GeodesicEquationSolver();
    // Test_Stencil(5);
    // Test_TriangleIntersection();
    // Test_SphericalBarycentricWeights();
    // Test_Quadrature(MyStencil(7));
    // Test_Quadrature(LebedevStencil(7));
    // Test_Quadrature(GaussLegendreStencil(7));
    // Test_SphericalHarmonicsExpansion();
    // Test_QuaternionRotation();
    // Test_IntensityAt();
    // Test_PhotonFourVelocity();
    // Test_MyAtan2();
    // Test_MySin_MyCos();
    // Test_Camera();
    // Test_Emission( 25, 45, 50,  20,    15,10);
    // Test_HarmonicsBenchmark();
    Test_UnstructuredMatrix();
    // Test_ConvexHull();
    // Test_GetTriangle();
    // Test_Radiation2(7);
    // Test_Radiation2(9);
    // Test_Radiation2(11);
    // Test_Radiation2(13);
    // Test_Radiation3SphereWave(9, StreamingType::FlatStatic);
    // Test_Radiation3SphereWave(9, StreamingType::FlatDynamic);
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