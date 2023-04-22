#include <iostream>
#include "../src/Stencil.h"
#include "../src/SphericalHarmonics.h"
#include "../src/Profiler.hh"
using namespace std;



void PrintStencilData(Stencil stencil)
{
    cout << "Stencil Data:" << endl;
    stencil.Print();
    cout << endl;
    cout << "Connected Triangles:" << endl;
    stencil.connectedTriangles.Print();
    cout << endl;
    cout << "Connected Vertices Order 1:" << endl;
    stencil.connectedVerticesOrder1.Print();
    cout << endl;
    cout << "Connected Vertices Order 2:" << endl;
    stencil.connectedVerticesOrder2.Print();
    cout << endl;
}



void Quadrature(const Stencil& stencil)
{
    cout << "nDir = " << stencil.nDir << endl;
    cout << "nOrder = " << stencil.nOrder << endl;
    cout << "nCoefficients = " << stencil.nCoefficients << endl;

    double* data = new double[stencil.nDir]();
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



void QuadratureRotated(const Stencil& stencil)
{
    cout << "nDir = " << stencil.nDir << endl;
    cout << "nOrder = " << stencil.nOrder << endl;
    cout << "nCoefficients = " << stencil.nCoefficients << endl;

    // Random direction:
	Tensor3 n;
    double norm = 2;
    while(norm > 1)
    {
        n[1] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
        n[2] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
        n[3] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
        norm = n.EuklNorm();
    }
    n = n.EuklNormalized();

    glm::vec3 from(0,0,1);
    glm::vec3 to(n[1],n[2],n[3]);
    glm::quat q(from,to);

    Tensor3 north(0,0,1);
    n.Print("n");
    north.Print("north");
    (q*north).Print("northRot");

    double data[stencil.nDir];
    for(size_t d=0; d<stencil.nDir; d++)
    {
        Tensor3 dir = q * stencil.C(d);
        data[d] = 0;
        for(int k=0; k<stencil.nCoefficients; k++)
            data[d] += (k+1) * SphericalHarmonicsXyz::Y(k,dir);
    }

    double coefficients[stencil.nDir];

    // GetCoefficients modified to use rotated stencil:
    for(size_t i=0; i<stencil.nCoefficients; i++)
        coefficients[i] = 0;
    for(size_t d=0; d<stencil.nDir; d++)
    {
        Tensor3 dir = q * stencil.C(d);
        double x = dir[1];
        double y = dir[2];
        double z = dir[3];
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * SphericalHarmonicsXyz::Y(i,x,y,z);
    }
    
    for(size_t i=0; i<stencil.nCoefficients; i++)
        cout << "c" << i << ": " << coefficients[i] << endl;
    cout << endl;

    cout << "Real Value, Compressed Value:" << endl;
    size_t d = stencil.nDir/2;
    cout << data[d] << ", " << SphericalHarmonicsXyz::GetValue(q*stencil.C(d),coefficients,stencil.nCoefficients) << endl << endl;
}



void QuaternionRotation()
{
    using namespace glm;
    {
        vec3 from(0,0,1);
        Tensor3 n = Tensor3(1,1,1).EuklNormalized();
        vec3 to(n[1],n[2],n[3]);
        quat q = quat(from, to);
        quat qInv = quat(to, from);

        GaussLegendreStencil stencil(15);

        ofstream file0(OUTPUTDIR + "QuaternionRotationNormal.txt");
        ofstream file1(OUTPUTDIR + "QuaternionRotationRotated.txt");
        file0 << "x, y, z, color" << endl;
        file1 << "x, y, z, color" << endl;
        for(size_t d=0; d<stencil.nDir; d++)
        {
            Tensor3 p = stencil.C(d);
            file0 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 0 << endl;
            p = q * p;
            file1 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 0 << endl;
        }
        file0.close();
        file1.close();
    }

    {
        vec3 from0(0,0,1);
        Tensor3 n0 = Tensor3(-1,0,1).EuklNormalized();
        vec3 to0(n0[1],n0[2],n0[3]);
        quat q0 = quat(from0, to0);
        
        vec3 from1(0,0,1);
        Tensor3 n1 = Tensor3(1,0,1).EuklNormalized();
        vec3 to1(n1[1],n1[2],n1[3]);
        quat q1 = quat(from1, to1);

        Tensor3 vL0(0.70710678,0,0.70710678);
        Tensor3 vW = q0 * vL0;
        Tensor3 vL1 = Invert(q1) * vW;

        vL0.Print("vL0");
        vW.Print("vW");
        vL1.Print("vL1");
        cout << endl;
    }
    
    cout << "2 files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



void HarmonicsBenchmark()
{
    LebedevStencil stencil(31);
    double* data1 = new double[stencil.nDir];
    double* data2 = new double[stencil.nDir];
    double* coefficients1 = new double[stencil.nCoefficients];
    double* coefficients2 = new double[stencil.nCoefficients];
    double errorTheshhold = 1e-4;

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
            // for(int k=0; k<stencil.nCoefficients; k++)
                // cout << "c[" << k << "]: " << coefficients1[k] << endl;
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
    cout << endl;
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
            // for(int k=0; k<stencil.nCoefficients; k++)
                // cout << "c[" << k << "]: " << coefficients2[k] << endl;
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



void ConvertStencilToUnityMesh(Stencil stencil, string name)
{
    // Setup vertices vector for convex hull triangulation:
    std::vector<Vector3> vertices;
    vertices.reserve(stencil.nDir);
    for(int d=0; d<stencil.nDir; d++)
    {
        Tensor3 c = stencil.C(d);
        Vector3 v(c[1],c[2],c[3]);
        vertices.push_back(v);
    }

    // Triangulate vertices:
    ConvexHull convexHull(vertices);
    convexHull.OriginalOrdering(vertices);

    convexHull.GetMesh().WriteToCsv(OUTPUTDIR, name + to_string(stencil.nOrder));
}


void NearestNeighbourGrid(Stencil stencil)
{
    ofstream file0(OUTPUTDIR + "nearestNeighbourGrid0.txt");
    ofstream file1(OUTPUTDIR + "nearestNeighbourGrid1.txt");
    ofstream file2(OUTPUTDIR + "nearestNeighbourGrid2.txt");
    file0 << "x, y, z, index\n";
    file1 << "x, y, z, index\n";
    file2 << "x, y, z, index\n";
    for(int j=0; j<stencil.sphereGrid.nPh; j++)
    for(int i=0; i<stencil.sphereGrid.nTh; i++)
    {
        Tensor3 c = stencil.sphereGrid.C(i,j);
        file0 << c[1] << ", " << c[2] << ", " << c[3] << ", " << stencil.GetNeighbourIndex0(c) << "\n";
        file1 << c[1] << ", " << c[2] << ", " << c[3] << ", " << stencil.GetNeighbourIndex1(c) << "\n";
        file2 << c[1] << ", " << c[2] << ", " << c[3] << ", " << stencil.GetNeighbourIndex2(c) << "\n";
    }
    file0.close();
    file1.close();
    file2.close();
    
    cout << "3 files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
}


void NearestNeighbourGridLoop(Stencil stencil)
{
    cout << "Startint infinite loop and testing random directions. Stop on fail. Cancel yourself." << endl;
    while(true)
    {
        // Random direction Vector:
        Tensor3 n;
        double norm = 2;
        while(norm > 1)
        {
            n[1] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
            n[2] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
            n[3] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
            norm = n.EuklNorm();
        }
        n = n.EuklNormalized();
        n.Print("n");

        // Nearest Vector via lookup:
        size_t d0 = stencil.GetNeighbourIndex0(n);
        size_t d1 = stencil.GetNeighbourIndex1(n);
        size_t d2 = stencil.GetNeighbourIndex2(n);
        cout << "d0: " << d0 << endl;
        cout << "d1: " << d1 << endl;
        cout << "d2: " << d2 << endl;
        double dot0 = Tensor3::Dot(n,stencil.C(d0));
        double dot1 = Tensor3::Dot(n,stencil.C(d1));
        double dot2 = Tensor3::Dot(n,stencil.C(d2));
        size_t d = -1;
        if      (dot0 >= dot1 && dot0 >= dot2) d = d0;
        else if (dot1 >= dot0 && dot1 >= dot2) d = d1;
        else    d = d2;
        cout << "d: " << d << endl;

        // Nearest Vector via brute force:
        size_t D;
        double maxDot = -1;
        for(size_t p=0; p<stencil.nDir; p++)
        {
            Tensor3 c = stencil.C(p);
            double dot = Tensor3::Dot(n,c);
            if(dot > maxDot)
            {
                maxDot = dot;
                D = p;
            }
        }
        cout << "D: " << D << endl;

        // Comparison:
        if(d != D)
        {
            size_t d0 = stencil.GetNeighbourIndex0(n);
            size_t d1 = stencil.GetNeighbourIndex1(n);
            size_t d2 = stencil.GetNeighbourIndex2(n);

            std::cout << " D: " << D  << std::endl;
            std::cout << " d: " << d  << std::endl;
            std::cout << "d0: " << d0 << std::endl;
            std::cout << "d1: " << d1 << std::endl;
            std::cout << "d2: " << d2 << std::endl;
            n.Print("n");
            stencil.C(D) .Print("c(D)");
            stencil.C(d) .Print("c(d)");
            stencil.C(d0).Print("c(d0)");
            stencil.C(d1).Print("c(d1)");
            stencil.C(d2).Print("c(d2)");
            std::cin.get();
        }
    }
}



int main()
{
    // PrintStencilData(GaussLegendreStencil(5));
    // PrintStencilData(LebedevStencil(9));
    // Quadrature(GaussLegendreStencil(7));
    // Quadrature(LebedevStencil(7));
    // Quadrature(GaussLegendreStencil(35));
    // Quadrature(LebedevStencil(47));
    // QuadratureRotated(LebedevStencil(7));
    // QuaternionRotation();
    // HarmonicsBenchmark();

    // for(int i=3; i<32; i+=2)
        // ConvertStencilToUnityMesh(LebedevStencil(i),"LebedevStencil");
    // ConvertStencilToUnityMesh(LebedevStencil(35),"LebedevStencil");
    // ConvertStencilToUnityMesh(LebedevStencil(41),"LebedevStencil");
    // ConvertStencilToUnityMesh(LebedevStencil(47),"LebedevStencil");
    // for(int i=3; i<36; i+=2)
        // ConvertStencilToUnityMesh(GaussLegendreStencil(i),"GaussLegendreStencil");

    // NearestNeighbourGrid(LebedevStencil(5));
    // NearestNeighbourGridLoop(LebedevStencil(11));
}