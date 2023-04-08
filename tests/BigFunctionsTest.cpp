#include <iostream>
#include "../src/Radiation.h"
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
    MyStencil stencil(15);
    LebedevStencil lebedevStencil(5);
    Camera camera;

    Radiation radiation(metric, stencil, lebedevStencil, camera, StreamingType::CurvedDynamic);
    radiation.UpdateSphericalHarmonicsCoefficients();
    
    ofstream file0((string)OUTPUTDIR + (string)"SphericalHarmonicsExpansionCoord.csv");
    ofstream file1((string)OUTPUTDIR + (string)"SphericalHarmonicsExpansionVeloc.csv");
    ofstream file2((string)OUTPUTDIR + (string)"GeodesicCoord.csv");
    ofstream file3((string)OUTPUTDIR + (string)"GeodesicVeloc.csv");
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

    cout << "4 files have been created in '" + (string)OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



void IntensityAt()
{
    // This test is simulating the Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF) function.
    MyStencil stencil(5);
    LebedevStencil lebedev(11);

    double I[stencil.nDir];
    double Inorth;
    double Isouth;

    ofstream file0((string)OUTPUTDIR + (string)"IntensityAtStencil.csv");
    ofstream file1((string)OUTPUTDIR + (string)"IntensityAtLebedev.csv");
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

    cout << "2 files have been created in '" + (string)OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



int main()
{
    SphericalHarmonicsExpansion();
    IntensityAt();
}