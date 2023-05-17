#include <iostream>
#include "../src/Radiation3.h"
using namespace std;



void SphericalHarmonicsExpansion()
{
    size_t nx, ny, nz;
    nx = ny = nz = 20;
    Coord start(0,0,0);
    Coord end(3,3,3);
    Grid grid(nx, ny, nz, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    GaussLegendreStencil stencil(15);
    LebedevStencil lebedevStencil(5);
    Camera camera;

    Radiation3 radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedDynamic);
    radiation.UpdateSphericalHarmonicsCoefficients();
    
    ofstream file0(OUTPUTDIR + "SphericalHarmonicsExpansionCoord.txt");
    ofstream file1(OUTPUTDIR + "SphericalHarmonicsExpansionVeloc.txt");
    ofstream file2(OUTPUTDIR + "GeodesicCoord.txt");
    ofstream file3(OUTPUTDIR + "GeodesicVeloc.txt");
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
        for(size_t d=0; d<stencil.nDir; d++)
        {
            // Get pos and vel by Spherical Harmonic Expansion:
            {
                Tensor3 direction = stencil.Ct3(d);
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
                Tensor3 c = stencil.Ct3(d);
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

    cout << "4 files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



double IntensityAtBarycentric(Tensor3 direction, const glm::quat& qSrc, Stencil& stencil, double* I)
{
    direction = Invert(qSrc) * direction;
    size_t d = stencil.NearestNeighbour(direction);
    Vector3Int triangle;
    Vector3 weights;
    Tensor3 rayOrigin(0,0,0);
    std::span triangles = stencil.TrianglesConnectedTo(d);
    for(auto t = begin(triangles); t != end(triangles); t++)
    {
        triangle = *t;
        Tensor3 v0 = stencil.Ct3(triangle[0]);
        Tensor3 v1 = stencil.Ct3(triangle[1]);
        Tensor3 v2 = stencil.Ct3(triangle[2]);
        if (BarycentricWeights(rayOrigin, direction, v0, v1, v2, weights))
            break;
    }
    double I0 = I[triangle[0]];
    double I1 = I[triangle[1]];
    double I2 = I[triangle[2]];
    return std::max(0.0, I0 * weights[0] + I1 * weights[1] + I2 * weights[2]);
}
double IntensityAtRayIntersection(Tensor3 direction, const glm::quat& qSrc, Stencil& stencil, double* I)
{
    direction = Invert(qSrc) * direction;
    size_t d = stencil.NearestNeighbour(direction);
    Vector3Int triangle;
    Tensor3 intersectionPoint;
    Tensor3 rayOrigin(0,0,0);
    std::span triangles = stencil.TrianglesConnectedTo(d);
    for(auto t = begin(triangles); t != end(triangles); t++)
    {
        triangle = *t;
        Tensor3 v0 = I[triangle[0]] * stencil.Ct3(triangle[0]);
        Tensor3 v1 = I[triangle[1]] * stencil.Ct3(triangle[1]);
        Tensor3 v2 = I[triangle[2]] * stencil.Ct3(triangle[2]);
        if (RayTriangleIntersection(rayOrigin, direction, v0, v1, v2, intersectionPoint))
            break;
    }
    return intersectionPoint.EuklNorm();
}
double IntensityAtSphereBarycentric(Tensor3 direction, const glm::quat& qSrc, Stencil& stencil, double* I)
{
    direction = Invert(qSrc) * direction;
    size_t d = stencil.NearestNeighbour(direction);
    Vector3Int triangle;
    Tensor3 weights;
    Tensor3 rayOrigin(0,0,0);
    std::span triangles = stencil.TrianglesConnectedTo(d);
    for(auto t = begin(triangles); t != end(triangles); t++)
    {
        triangle = *t;
        Tensor3 v0 = stencil.Ct3(triangle[0]);
        Tensor3 v1 = stencil.Ct3(triangle[1]);
        Tensor3 v2 = stencil.Ct3(triangle[2]);
        if (SphericalBarycentricWeights(direction, v0, v1, v2, weights))
            break;
    }
    double I0 = I[triangle[0]];
    double I1 = I[triangle[1]];
    double I2 = I[triangle[2]];
    return std::max(0.0, I0 * weights[1] + I1 * weights[2] + I2 * weights[3]);
}
double IntensityAtInvDist(Tensor3 direction, const glm::quat& qSrc, Stencil& stencil, double* I)
{
    return 0;
    //direction = Invert(qSrc) * direction;
    //size_t d = stencil.GetNearestNeighbour(direction);
    //double value = 0;
    //double invDistSum = 0;
    //for(size_t p=stencil.connectedVerticesOrder1.Start(d); p<stencil.connectedVerticesOrder1.End(d); p++)
    //{
    //    size_t index = stencil.connectedVerticesOrder1[p];
    //    double dist = Tensor3::UnitSphereNorm(direction,stencil.C(index));
    //    if(dist < 1e-16)    // vTempIF == stencil.C(index)
    //        return I[index];
    //    double invDist = 1.0 / dist;
    //    value += I[index] * invDist;
    //    invDistSum += invDist;
    //}
    //return std::max(0.0, value/invDistSum);
}
void SphereInterpolation(Stencil srcStencil, Stencil dstStencil)
{
    // Rotation of source stencil:
    Tensor3 fSrc = Tensor3(0,0,1).EuklNormalized();
    Tensor3 tSrc   = Tensor3(1,1,1).EuklNormalized();
    glm::vec3 fromSrc(fSrc[1],fSrc[2],fSrc[3]);
    glm::vec3 toSrc  (tSrc[1],tSrc[2],tSrc[3]);
    glm::quat qSrc(fromSrc,toSrc);
    
    // Fill intensities of source stencil:
    double I0[srcStencil.nDir];
    double sigma = 2;
    for(int d=0; d<srcStencil.nDir; d++)
    {
        Tensor3 c = qSrc * srcStencil.Ct3(d);
        I0[d] = exp(sigma * Tensor3::Dot(tSrc, c));
    }

    // Write source stencil to file:
    ofstream file0(OUTPUTDIR + "srcStencil.txt");
    file0 << "x,y,z,color\n";
    for(int d=0; d<srcStencil.nDir; d++)
    {
        Tensor3 c = I0[d] * (qSrc * srcStencil.Ct3(d));
        file0 << c[1] << "," << c[2] << "," << c[3] << "," << 0 << "\n";
    }
    file0.close();



    // Rotation of destination stencil:
    Tensor3 fDst = Tensor3(0,0,1).EuklNormalized();
    Tensor3 tDst   = Tensor3(-1,-1,-1).EuklNormalized();
    glm::vec3 fromDst(fDst[1],fDst[2],fDst[3]);
    glm::vec3 toDst  (tDst[1],tDst[2],tDst[3]);
    glm::quat qDst(fromDst,toDst);

    // Fill intensities of destination stencil:
    double I1[dstStencil.nDir];
    double I2[dstStencil.nDir];
    double I3[dstStencil.nDir];
    double I4[dstStencil.nDir];
    for(int d=0; d<dstStencil.nDir; d++)
    {
        Tensor3 c = qDst * dstStencil.Ct3(d);
        I1[d] = IntensityAtBarycentric(c, qSrc, srcStencil, I0);
        I2[d] = IntensityAtRayIntersection(c, qSrc, srcStencil, I0);
        I3[d] = IntensityAtSphereBarycentric(c, qSrc, srcStencil, I0);
        I4[d] = IntensityAtInvDist(c, qSrc, srcStencil, I0);
    }

    // Write destination stencil to file:
    ofstream file1(OUTPUTDIR + "dstStencilBarycentric.txt");
    ofstream file2(OUTPUTDIR + "dstStencilSphereRayIntersection.txt");
    ofstream file3(OUTPUTDIR + "dstStencilSphereBarycentric.txt");
    ofstream file4(OUTPUTDIR + "dstStencilInvDist.txt");
    file1 << "x,y,z,color\n";
    file2 << "x,y,z,color\n";
    file3 << "x,y,z,color\n";
    file4 << "x,y,z,color\n";
    for(int d=0; d<dstStencil.nDir; d++)
    {
        Tensor3 c1 = I1[d] * (qDst * dstStencil.Ct3(d));
        Tensor3 c2 = I2[d] * (qDst * dstStencil.Ct3(d));
        Tensor3 c3 = I3[d] * (qDst * dstStencil.Ct3(d));
        Tensor3 c4 = I4[d] * (qDst * dstStencil.Ct3(d));
        file1 << c1[1] << "," << c1[2] << "," << c1[3] << "," << 0.25 << "\n";
        file2 << c2[1] << "," << c2[2] << "," << c2[3] << "," << 0.50 << "\n";
        file3 << c3[1] << "," << c3[2] << "," << c3[3] << "," << 0.75 << "\n";
        file4 << c4[1] << "," << c4[2] << "," << c4[3] << "," << 1.00 << "\n";
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();

    cout << "5 files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



int main()
{
    // SphericalHarmonicsExpansion();
    // SphereInterpolation(GaussLegendreStencil(15));
    SphereInterpolation(LebedevStencil(15), GaussLegendreStencil(21));
    // SphereInterpolation(LebedevStencil(35), GaussLegendreStencil(21));
    // SphereInterpolation(GaussLegendreStencil(5), GaussLegendreStencil(15));
    // SphereInterpolation(LebedevStencil(11));
}