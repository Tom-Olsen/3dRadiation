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
    cout << "Connected Vertices:" << endl;
    stencil.connectedVertices.Print();
    cout << endl;
}



void Quadrature(const Stencil& stencil)
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
    vec3 from(0,0,1);
    vec3 to(1,0,0);
    quat q = quat(from, to);
    quat qInv = quat(to, from);

    MyStencil stencil(15);

    ofstream file0((string)OUTPUTDIR + (string)"QuaternionRotationNormal.csv");
    ofstream file1((string)OUTPUTDIR + (string)"QuaternionRotationRotated.csv");
    file0 << "x, y, z, color" << endl;
    file1 << "x, y, z, color" << endl;
    
    for(size_t d1=0; d1<stencil.nPh; d1++)
    for(size_t d0=0; d0<stencil.nTh; d0++)
    {// Bulk:
        Tensor3 p = stencil.C(d0,d1);
        file0 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 0 << endl;
        p = q * p;
        file1 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 0 << endl;
    }
    {// North Pole:
        Tensor3 p(0,0,1);
        file0 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 1 << endl;
        p = q * p;
        file1 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 1 << endl;
    }

    file0.close();
    file1.close();

    Tensor3 p(1,2,3);
    Tensor3 pRot = q * p;
    Tensor3 pRotInvRot = qInv * pRot;
    Tensor3 pRotInvRot2 = Invert(q) * pRot;

    qPrint(q,"q");
    qPrint(qInv,"qInv0");
    qPrint(Invert(q),"qInv1");
    cout << endl;
    p.Print("   p");
    pRot.Print("pRot");
    pRotInvRot.Print("  p?");
    pRotInvRot2.Print("  p?");
    
    cout << "2 files have been created in '" + (string)OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
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



int main()
{
    PrintStencilData(LebedevStencil(9));
    // Quadrature(MyStencil(7));
    // Quadrature(GaussLegendreStencil(7));
    // Quadrature(LebedevStencil(7));
    // QuadratureRotated(LebedevStencil(7));
    // QuaternionRotation();
    // HarmonicsBenchmark();
}