#include <iostream>
#include "../src/Stencil.h"
#include "../src/SphericalHarmonics.h"
#include "../src/Profiler.hh"
using namespace std;



void PrintStencilData(Stencil stencil)
{
    // Test Point:
    double theta = 1.37;
    double phi = 0.86;
    Tensor3 p(MySin(theta)*MyCos(phi), MySin(theta)*MySin(phi), MyCos(theta));
    // Tensor3 p(1,0,0);

    cout << "Stencil Data:" << endl;
    stencil.Print();
    cout << endl;

    cout << "Connected Triangles:" << endl;
    for(int d=0; d<stencil.nDir; d++)
        PrintList(stencil.TrianglesConnectedTo(d), to_string(d));
    cout << endl;

    cout << "Voronoi Cells:" << endl;
    for(int d=0; d<stencil.nDir; d++)
        PrintList(stencil.VoronoiCellOf(d), to_string(d));
    cout << endl;

    cout << "Voronoi Neighbours:" << endl;
    for(int d=0; d<stencil.nDir; d++)
        PrintList(stencil.VoronoiNeighbourOf(d), to_string(d));
    cout << endl;
    
    // Change access modifier of Stencil::VirtualVoronoiCellOf to test this!
    // double theta = 1.37;
    // double phi = 0.86;
    // Tensor3 p(MySin(theta)*MyCos(phi), MySin(theta)*MySin(phi), MyCos(theta));
    // std::vector<size_t> naturalNeighbours;
    // std::vector<Vector3> pointsInsideCell;
    // cout << "Virtual Voronoi Cell:" << endl;
    // cout << "p = " << p << endl;
    // std::vector<Vector3> virtualCell = stencil.VirtualVoronoiCellOf(p, naturalNeighbours, pointsInsideCell);
    // PrintList(virtualCell, "Virtual Cell");
    // cout << endl;

    cout << "Voronoi Weights and Supports:" << endl;
    std::vector<size_t> naturalNeighbours;
    std::vector<double> voronoiWeights = stencil.VoronoiWeights(p, naturalNeighbours);
    PrintList(voronoiWeights, "weights");
    PrintList(naturalNeighbours, "naturalNeighbours");
    cout << endl;

    cout << "Barycentric Weights and Support:" << endl;
    Vector3Int triangle;
    Vector3 BaryWeights = stencil.BarycentricWeights(p,triangle);
    cout << "weights: " << BaryWeights << endl;
    cout << "triangle: " << triangle << endl;
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
        Tensor3 dir = stencil.Ct3(d);
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
    cout << data[d] << ", " << SphericalHarmonicsXyz::GetValue(stencil.Ct3(d),coefficients,stencil.nCoefficients) << endl << endl;
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
        Tensor3 dir = q * stencil.Ct3(d);
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
        Tensor3 dir = q * stencil.Ct3(d);
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
    cout << data[d] << ", " << SphericalHarmonicsXyz::GetValue(q*stencil.Ct3(d),coefficients,stencil.nCoefficients) << endl << endl;
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



void ConvertStencilToUnityMesh(Stencil stencil)
{
    // Setup vertices vector for convex hull triangulation:
    std::vector<Vector3> vertices;
    vertices.reserve(stencil.nDir);
    for(int d=0; d<stencil.nDir; d++)
    {
        Tensor3 c = stencil.Ct3(d);
        Vector3 v(c[1],c[2],c[3]);
        vertices.push_back(v);
    }

    // Triangulate vertices:
    ConvexHull convexHull(vertices);
    convexHull.OriginalOrdering(vertices);

    convexHull.GetMesh().WriteToCsv(OUTPUTDIR, stencil.name);
}



int main()
{
    // PrintStencilData(GaussLegendreStencil(5));
    // PrintStencilData(LebedevStencil(3));
    PrintStencilData(LebedevStencil(3,2,4,M_PI/8.0));
    // Quadrature(GaussLegendreStencil(7));
    // Quadrature(LebedevStencil(7));
    // Quadrature(GaussLegendreStencil(35));
    // Quadrature(LebedevStencil(47));
    // QuadratureRotated(LebedevStencil(7));
    // HarmonicsBenchmark();

    // for(int i=3; i<32; i+=2)
        // ConvertStencilToUnityMesh(LebedevStencil(i));
    // ConvertStencilToUnityMesh(LebedevStencil(35));
    // ConvertStencilToUnityMesh(LebedevStencil(41));
    // ConvertStencilToUnityMesh(LebedevStencil(47));
    // for(int i=3; i<36; i+=2)
        // ConvertStencilToUnityMesh(GaussLegendreStencil(i));
}