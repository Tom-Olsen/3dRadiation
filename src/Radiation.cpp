#include "Radiation.h"

std::string StreamingName(int n)
{
	std::string s("unknown");
	switch (n)
	{
   		case 0: { s = "FlatStatic";	     } break;
   		case 1: { s = "FlatDynamic";	 } break;
   		case 2: { s = "GeodesicStatic";	 } break;
   		case 3: { s = "GeodesicDynamic"; } break;
   		case 4: { s = "CurvedStatic";    } break;
   		case 5: { s = "CurvedDynamic";	 } break;
		default: { exit_on_error("Invalid StreamingType"); }
	}
	return s;
}



Radiation::Radiation(Metric& metric_, Stencil& stencil_, LebedevStencil& lebedevStencil_, StreamingType streamingType_):
grid(metric_.grid), metric(metric_), stencil(stencil_), lebedevStencil(lebedevStencil_), streamingType(streamingType_)
{
	isInitialGridPoint = new bool[grid.nxyz]();
	initialE           = new double[grid.nxyz]();
	initialNx          = new double[grid.nxyz]();
	initialNy          = new double[grid.nxyz]();
	initialNz          = new double[grid.nxyz]();
	initialKappa0      = new double[grid.nxyz]();
	initialKappa1      = new double[grid.nxyz]();
	initialKappaA      = new double[grid.nxyz]();
	initialEta         = new double[grid.nxyz]();
	nx     = new double[grid.nxyz]();
	ny     = new double[grid.nxyz]();
	nz     = new double[grid.nxyz]();
	nxNew  = new double[grid.nxyz]();
	nyNew  = new double[grid.nxyz]();
	nzNew  = new double[grid.nxyz]();
	E      = new double[grid.nxyz]();
	Fx     = new double[grid.nxyz]();
	Fy     = new double[grid.nxyz]();
	Fz     = new double[grid.nxyz]();
	Pxx    = new double[grid.nxyz]();
	Pxy    = new double[grid.nxyz]();
	Pxz    = new double[grid.nxyz]();
	Pyy    = new double[grid.nxyz]();
	Pyz    = new double[grid.nxyz]();
	Pzz    = new double[grid.nxyz]();
	E_LF   = new double[grid.nxyz]();
	Fx_LF  = new double[grid.nxyz]();
	Fy_LF  = new double[grid.nxyz]();
	Fz_LF  = new double[grid.nxyz]();
	Pxx_LF = new double[grid.nxyz]();
	Pxy_LF = new double[grid.nxyz]();
	Pxz_LF = new double[grid.nxyz]();
	Pyy_LF = new double[grid.nxyz]();
	Pyz_LF = new double[grid.nxyz]();
	Pzz_LF = new double[grid.nxyz]();
	kappa0 = new double[grid.nxyz]();
	kappa1 = new double[grid.nxyz]();
	kappaA = new double[grid.nxyz]();
	eta    = new double[grid.nxyz]();
	I      = new double[grid.nxyz * stencil.nThPh]();
	Inew   = new double[grid.nxyz * stencil.nThPh]();
	Inorth = new double[grid.nxyz]();
	Isouth = new double[grid.nxyz]();
	coefficientsS  = new double[grid.nxyz * lebedevStencil.nDir];
	coefficientsX  = new double[grid.nxyz * lebedevStencil.nDir];
	coefficientsY  = new double[grid.nxyz * lebedevStencil.nDir];
	coefficientsZ  = new double[grid.nxyz * lebedevStencil.nDir];
	coefficientsCx = new double[grid.nxyz * lebedevStencil.nDir];
	coefficientsCy = new double[grid.nxyz * lebedevStencil.nDir];
	coefficientsCz = new double[grid.nxyz * lebedevStencil.nDir];
}



Radiation::~Radiation()
{
	delete[] isInitialGridPoint;
	delete[] initialE;
	delete[] initialNx;
	delete[] initialNy;
	delete[] initialNz;
	delete[] initialKappa0;
	delete[] initialKappa1;
	delete[] initialKappaA;
	delete[] initialEta;
	delete[] nx;
	delete[] ny;
	delete[] nz;
	delete[] nxNew;
	delete[] nyNew;
	delete[] nzNew;
	delete[] E;
	delete[] Fx;
	delete[] Fy;
	delete[] Fz;
	delete[] Pxx;
	delete[] Pxy;
	delete[] Pxz;
	delete[] Pyy;
	delete[] Pyz;
	delete[] Pzz;
	delete[] E_LF;
	delete[] Fx_LF;
	delete[] Fy_LF;
	delete[] Fz_LF;
	delete[] Pxx_LF;
	delete[] Pxy_LF;
	delete[] Pxz_LF;
	delete[] Pyy_LF;
	delete[] Pyz_LF;
	delete[] Pzz_LF;
	delete[] kappa0;
	delete[] kappa1;
	delete[] kappaA;
	delete[] eta;
	delete[] I;
	delete[] Inew;
	delete[] Inorth;
	delete[] Isouth;
	delete[] coefficientsS;
	delete[] coefficientsX;
	delete[] coefficientsY;
	delete[] coefficientsZ;
	delete[] coefficientsCx;
	delete[] coefficientsCy;
	delete[] coefficientsCz;
}



void Radiation::NormalizeInitialDirections()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		Tensor3 n(initialNx[ijk], initialNy[ijk], initialNz[ijk]);
		double norm = n.EuklNorm();

		if(norm < normThreshhold)
		{
			initialNx[ijk] = 0;
			initialNy[ijk] = 0;
			initialNz[ijk] = 0;
		}
		else
		{
			initialNx[ijk] = n[1] / norm;
			initialNy[ijk] = n[2] / norm;
			initialNz[ijk] = n[3] / norm;
		}
	}
}



void Radiation::LoadInitialData()
{
	PROFILE_FUNCTION();
	bool isDynamicStreaming = (streamingType == StreamingType::FlatDynamic || streamingType == StreamingType::GeodesicDynamic || streamingType == StreamingType::CurvedDynamic);
	srand((unsigned) time(NULL));

	#pragma omp parallel for
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		if(isInitialGridPoint[ijk])
		{
			nx[ijk] = initialNx[ijk];
			ny[ijk] = initialNy[ijk];
			nz[ijk] = initialNz[ijk];
			kappa0[ijk] = initialKappa0[ijk];
			kappa1[ijk] = initialKappa1[ijk];
			kappaA[ijk] = initialKappaA[ijk];
			eta[ijk] = initialEta[ijk];

			if(nx[ijk] == 0 && ny[ijk] == 0 && nz[ijk] == 0)
			{// Uniform intensity distribution:
				if (isDynamicStreaming)
				{// Random direction for uniformity:
					Tensor3 n(0.0);
					n[1] = 2.0 * ((float) rand() / RAND_MAX) - 1.0;
					n[2] = 2.0 * ((float) rand() / RAND_MAX) - 1.0;
					n[3] = 2.0 * ((float) rand() / RAND_MAX) - 1.0;
					n = n.EuklNormalized();
					nx[ijk] = n[1];
					ny[ijk] = n[2];
					nz[ijk] = n[3];
				}
				else
				{// (0,0,0) is an invalid direction. Default direction is towards z:
					nx[ijk] = 0;	
					ny[ijk] = 0;	
					nz[ijk] = 1;	
				}
				for(int d=0; d<stencil.nThPh; d++)
					I[ijk + d*grid.nxyz] = initialE[ijk];
			}
			else
			{// Kent intensity distribution:
			 // https://en.wikipedia.org/wiki/Kent_distribution
			 	Tensor3 n = (isDynamicStreaming) ? Tensor3(0,0,1) : Tensor3(nx[ijk],ny[ijk],nz[ijk]);
				
				// In case of dynamic streaming the stencil is rotated such that its north pole points towards n, whenever cxyz is used.
				// This means that the Kent distribution needs to point north for the dynamic streaming case.

				for(int d1=0; d1<stencil.nPh; d1++)
				for(int d0=0; d0<stencil.nTh; d0++)
				{
					int d = stencil.Index(d0,d1);
					Tensor3 p = stencil.Cxyz(d0,d1);
					I[ijk + d*grid.nxyz] = initialE[ijk] * exp(stencil.sigma * Tensor3::Dot(n, p)) / (4.0 * M_PI * sinh(stencil.sigma));
				}
			}
		}
	}
}



void Radiation::NormalizeInitialData()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int ijk=0; ijk<grid.nxyz; ijk++)
		if(isInitialGridPoint[ijk])
			initialE[ijk] *= initialE[ijk] / E_LF[ijk];
}



void Radiation::UpdateSphericalHarmonicsCoefficients()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
	{
		int ijk = grid.Index(i,j,k);
		double dataS[lebedevStencil.nDir];
		double dataX[lebedevStencil.nDir];
		double dataY[lebedevStencil.nDir];
		double dataZ[lebedevStencil.nDir];
    	double dataCx[lebedevStencil.nDir];
    	double dataCy[lebedevStencil.nDir];
    	double dataCz[lebedevStencil.nDir];
		Coord xyz0 = grid.xyz(i,j,k);
		double alpha = metric.GetAlpha(ijk);

		for(int d=0; d<lebedevStencil.nDir; d++)
    	{
			// Initial data for geodesic equation:
			double s = 1;
			Coord xyz = xyz0;
            Tensor3 c = lebedevStencil.Cxyz(d);
            Tensor4 u(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
            Tensor3 v = Vec3ObservedByEulObs<IF,LF>(u, xyz, metric);

			// Solve geodesic equation backwards:
			if(!metric.InsideBH(xyz))
			{
    	    	// s *= Euler_GeodesicEquation<-1>(grid.dt,x,v,metric);
    	    	s *= RK45_GeodesicEquation(grid.dt, xyz, v, metric);
			}
			else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
				v = Tensor3(0.0);

			// Final data points for fourier expansion:
			dataS[d] = 1.0/s;
			dataX[d] = xyz[1];
			dataY[d] = xyz[2];
			dataZ[d] = xyz[3];
			dataCx[d] = v[1];
			dataCy[d] = v[2];
			dataCz[d] = v[3];
    	}
		std::vector<double> cS  = SphericalHarmonics::GetCoefficients(lebedevStencil, dataS , lebedevStencil.nCoefficients);
		std::vector<double> cX  = SphericalHarmonics::GetCoefficients(lebedevStencil, dataX , lebedevStencil.nCoefficients);
		std::vector<double> cY  = SphericalHarmonics::GetCoefficients(lebedevStencil, dataY , lebedevStencil.nCoefficients);
		std::vector<double> cZ  = SphericalHarmonics::GetCoefficients(lebedevStencil, dataZ , lebedevStencil.nCoefficients);
		std::vector<double> cCx = SphericalHarmonics::GetCoefficients(lebedevStencil, dataCx, lebedevStencil.nCoefficients);
		std::vector<double> cCy = SphericalHarmonics::GetCoefficients(lebedevStencil, dataCy, lebedevStencil.nCoefficients);
		std::vector<double> cCz = SphericalHarmonics::GetCoefficients(lebedevStencil, dataCz, lebedevStencil.nCoefficients);

		for(int f=0; f<lebedevStencil.nCoefficients; f++)
		{
			int fijk = f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy;
			coefficientsS[fijk] = cS[f];
			coefficientsX[fijk] = cX[f];
			coefficientsY[fijk] = cY[f];
			coefficientsZ[fijk] = cZ[f];
			coefficientsCx[fijk] = cCx[f];
			coefficientsCy[fijk] = cCy[f];
			coefficientsCz[fijk] = cCz[f];
		}
	}
}



void Radiation::ComputeMomentsIF()
{
	PROFILE_FUNCTION();
	constexpr double fourPiInv = 1.0 / (4.0 * M_PI);
	
	#pragma omp parallel for
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		E[ijk]   = 0.0;
		Fx[ijk]  = 0.0;
		Fy[ijk]  = 0.0;
		Fy[ijk]  = 0.0;
		Pxx[ijk] = 0.0;
		Pxy[ijk] = 0.0;
		Pxz[ijk] = 0.0;
		Pyy[ijk] = 0.0;
		Pyz[ijk] = 0.0;
		Pzz[ijk] = 0.0;
		glm::vec3 from(0,0,1);
		glm::vec3 to(nx[ijk],ny[ijk],nz[ijk]);
		glm::quat q(from,to);
		for(int d1=0; d1<stencil.nPh; d1++)
		for(int d0=0; d0<stencil.nTh; d0++)
		{
			Tensor3 cxyz = q * stencil.Cxyz(d0,d1);
			int d = stencil.Index(d0,d1);
			E[ijk]   += stencil.W(d0,d1) * I[ijk + d*grid.nxyz];
			Fx[ijk]  += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[1];
			Fy[ijk]  += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[2];
			Fz[ijk]  += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[3];
			Pxx[ijk] += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[1] * cxyz[1];
			Pxy[ijk] += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[1] * cxyz[2];
			Pxz[ijk] += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[1] * cxyz[3];
			Pyy[ijk] += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[2] * cxyz[2];
			Pyz[ijk] += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[2] * cxyz[3];
			Pzz[ijk] += stencil.W(d0,d1) * I[ijk + d*grid.nxyz] * cxyz[3] * cxyz[3];
		}
        E[ijk]   *= fourPiInv;
        Fx[ijk]  *= fourPiInv;
        Fy[ijk]  *= fourPiInv;
        Fz[ijk]  *= fourPiInv;
        Pxx[ijk] *= fourPiInv;
        Pxy[ijk] *= fourPiInv;
        Pxz[ijk] *= fourPiInv;
        Pyy[ijk] *= fourPiInv;
        Pyz[ijk] *= fourPiInv;
        Pzz[ijk] *= fourPiInv;
	}
}
void Radiation::ComputeMomentsLF()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		Tensor4x4 EnergyMomentumTensorIF
		( E[ijk], Fx[ijk], Fy[ijk], Fz[ijk],
		 Fx[ijk],Pxx[ijk],Pxy[ijk],Pxz[ijk],
		 Fy[ijk],Pxy[ijk],Pyy[ijk],Pyz[ijk],
		 Fz[ijk],Pxz[ijk],Pyz[ijk],Pzz[ijk]);
		Tensor4x4 EnergyMomentumTensorLF = TransformIFtoLF(EnergyMomentumTensorIF, metric.GetTetrad(ijk));

		E_LF[ijk]   = EnergyMomentumTensorLF[{0,0}];
		Fx_LF[ijk]  = EnergyMomentumTensorLF[{0,1}];
		Fy_LF[ijk]  = EnergyMomentumTensorLF[{0,2}];
		Fz_LF[ijk]  = EnergyMomentumTensorLF[{0,3}];
		Pxx_LF[ijk] = EnergyMomentumTensorLF[{1,1}];
		Pxy_LF[ijk] = EnergyMomentumTensorLF[{1,2}];
		Pxz_LF[ijk] = EnergyMomentumTensorLF[{1,3}];
		Pyy_LF[ijk] = EnergyMomentumTensorLF[{2,2}];
		Pyz_LF[ijk] = EnergyMomentumTensorLF[{2,3}];
		Pzz_LF[ijk] = EnergyMomentumTensorLF[{3,3}];
	}
}



void Radiation::SetIntensitiesNorthSouth()
{
	PROFILE_FUNCTION();
	int d0north = 0;
	int d0south = stencil.nTh - 1;
	#pragma omp parallel for
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		Inorth[ijk] = Isouth[ijk] = 0;
		for(int d1=0; d1<stencil.nPh; d1++)
		{
			Inorth[ijk] += I[ijk + stencil.Index(d0north,d1) * grid.nxyz];
			Isouth[ijk] += I[ijk + stencil.Index(d0south,d1) * grid.nxyz];
		}
		Inorth[ijk] /= stencil.nPh;
		Isouth[ijk] /= stencil.nPh;
	}
	//for(int k = 0; k < metric.grid.nz; k++)
	//{
	//	for(int j = 0; j < metric.grid.ny; j++)
	//	{
	//		for(int i = 0; i < metric.grid.nx; i++)
	//		{
	//			int ijk = grid.Index(i,j,k);
	//			std::cout << Format(Inorth[ijk],2) << ", ";
	//		}
	//		std::cout << std::endl;
	//	}
	//	std::cout << std::endl;
	//	std::cout << std::endl;
	//}
}



Coord Radiation::GetTempCoordinate(int i, int j, int k, double theta, double phi)
{
	std::vector<double> cx(lebedevStencil.nCoefficients);
	std::vector<double> cy(lebedevStencil.nCoefficients);
	std::vector<double> cz(lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
	{
		cx[f] = coefficientsX[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
		cy[f] = coefficientsY[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
		cz[f] = coefficientsZ[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
	}
	return Coord(SphericalHarmonics::GetValue(theta,phi,cx,lebedevStencil.nCoefficients),
				 SphericalHarmonics::GetValue(theta,phi,cy,lebedevStencil.nCoefficients),
				 SphericalHarmonics::GetValue(theta,phi,cz,lebedevStencil.nCoefficients));
}
Tensor3 Radiation::GetTemp3Velocity(int i, int j, int k, double theta, double phi)
{
	std::vector<double> cCx(lebedevStencil.nCoefficients);
	std::vector<double> cCy(lebedevStencil.nCoefficients);
	std::vector<double> cCz(lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
	{
		cCx[f] = coefficientsCx[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
		cCy[f] = coefficientsCy[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
		cCz[f] = coefficientsCz[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
	}
	return Tensor3(SphericalHarmonics::GetValue(theta,phi,cCx,lebedevStencil.nCoefficients),
				   SphericalHarmonics::GetValue(theta,phi,cCy,lebedevStencil.nCoefficients),
				   SphericalHarmonics::GetValue(theta,phi,cCz,lebedevStencil.nCoefficients));
}
double Radiation::GetFrequencyShift(int i, int j, int k, double theta, double phi)
{
	std::vector<double> cs(lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		cs[f] = coefficientsS[f + i*lebedevStencil.nCoefficients + j*lebedevStencil.nCoefficients*grid.nx + k*lebedevStencil.nCoefficients*grid.nxy];
	return SphericalHarmonics::GetValue(theta,phi,cs,lebedevStencil.nCoefficients);
}



double Radiation::IntensityAt(int ijk, Tensor3 vTempIF)
{
	// vTempIF is given in world space. Transform it to local space by applying inverse quaternion:
	glm::vec3 from(0,0,1);
	glm::vec3 to(nx[ijk],ny[ijk],nz[ijk]);
	glm::quat q(to,from);	// from and to are swapped to get inverse rotation.
	vTempIF = q * vTempIF;

	// Index on (old) sphere grid of point (theta,phi).
	double th = stencil.d0(vTempIF.Theta());
	double ph = stencil.d1(vTempIF.Phi());

	// Indices of nearest directions on sphere. Note that phi is cyclic and theta is not.
    int th0 = floor(th);
    int th1 = th0 + 1;
    int ph0 = ((int)floor(ph) + stencil.nPh) % stencil.nPh;
    int ph1 = (ph0 + 1) % stencil.nPh;

	// Get intensities at nearest directions. Use north pole for t0==-1 and south pole for t1==nTh.
    double I00 = (th0 == -1)          ? Inorth[ijk] : I[ijk + stencil.Index(th0,ph0) * grid.nxyz];
    double I01 = (th0 == -1)          ? Inorth[ijk] : I[ijk + stencil.Index(th0,ph1) * grid.nxyz];
    double I10 = (th1 == stencil.nTh) ? Isouth[ijk] : I[ijk + stencil.Index(th1,ph0) * grid.nxyz];
    double I11 = (th1 == stencil.nTh) ? Isouth[ijk] : I[ijk + stencil.Index(th1,ph1) * grid.nxyz];

	// Fractional part of sphere grid index. Near poles it must be streched:
	// North: from [0.5, 1.0] to [0.0, 1.0]
	// South: from [0.0, 0.5] to [0.0, 1.0]
	double tfrac = th - floor(th);
	double pfrac = ph - floor(ph);
    if (th0 == -1)			// North
        tfrac = 2.0 * (tfrac - 0.5);
    if (th1 == stencil.nTh)	// South
        tfrac *= 2.0;

	double value = BilinearInterpolation(tfrac, pfrac, I00, I01, I10, I11);
	return std::max(value,0.0);
}



Tensor3 Radiation::AverageF(int i, int j, int k)
{
	Tensor3 averageF(0.0);
	for(int c=-1; c<=1; c++)
	for(int b=-1; b<=1; b++)
	for(int a=-1; a<=1; a++)
	{
		int index = grid.Index(i+a, j+b, k+c);
		averageF[1] += Fx[index];
		averageF[2] += Fy[index];
		averageF[3] += Fz[index];
	}
	averageF[1] /= 27.0;
	averageF[2] /= 27.0;
	averageF[3] /= 27.0;
	return averageF;
}



void Radiation::StreamFlatStatic()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int d1=0; d1<stencil.nPh; d1++)
	for(int d0=0; d0<stencil.nTh; d0++)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
	{
		int d = stencil.Index(d0,d1);	// Index of direction d
		int ijk = grid.Index(i,j,k);	// Index of lattice point ijk
		int ijkd = grid.Index(i,j,k,d);	// Index of population d at lattice point ijk
		
		// Flat Streaming:
		Tensor3 vTempIF = stencil.Cxyz(d0,d1);
		Coord xyzTemp = grid.xyz(i,j,k);
		xyzTemp[1] -= vTempIF[1] * grid.dt;
		xyzTemp[2] -= vTempIF[2] * grid.dt;
		xyzTemp[3] -= vTempIF[3] * grid.dt;

		// Get 8 nearest Grid Points:
		double iTemp = grid.i(xyzTemp[1]);
		double jTemp = grid.j(xyzTemp[2]);
		double kTemp = grid.k(xyzTemp[3]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;
		int k0 = std::floor(kTemp);	int k1 = k0 + 1;

		// Interpolate intensity from neighbouring 8 lattice points to temporary point:
		Inew[ijkd]
		= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 I[grid.Index(i0,j0,k0,d)], I[grid.Index(i0,j0,k1,d)], I[grid.Index(i0,j1,k0,d)], I[grid.Index(i0,j1,k1,d)],
		 I[grid.Index(i1,j0,k0,d)], I[grid.Index(i1,j0,k1,d)], I[grid.Index(i1,j1,k0,d)], I[grid.Index(i1,j1,k1,d)]);
	}
	std::swap(I,Inew);
}
void Radiation::StreamFlatDynamic()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int d1=0; d1<stencil.nPh; d1++)
	for(int d0=0; d0<stencil.nTh; d0++)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
	{
		int d = stencil.Index(d0,d1);	// Index of direction d
		int ijk = grid.Index(i,j,k);	// Index of lattice point ijk
		int ijkd = grid.Index(i,j,k,d);	// Index of population d at lattice point ijk
		
		// Find new stencil direction:
		Tensor3 averageF = AverageF(i,j,k);
		Tensor3 n = (averageF.EuklNorm()>normThreshhold) ? averageF.EuklNormalized() : Tensor3(nx[ijk],ny[ijk],nz[ijk]);
		nxNew[ijk] = n[1];
		nyNew[ijk] = n[2];
		nzNew[ijk] = n[3];

		// Flat Streaming:
		glm::vec3 from(0,0,1);
		glm::vec3 to(nxNew[ijk],nyNew[ijk],nzNew[ijk]);
		glm::quat q(from,to);
		Tensor3 vTempIF = q * stencil.Cxyz(d0,d1);
		Coord xyzTemp = grid.xyz(i,j,k);
		xyzTemp[1] -= vTempIF[1] * grid.dt;
		xyzTemp[2] -= vTempIF[2] * grid.dt;
		xyzTemp[3] -= vTempIF[3] * grid.dt;

		// Get 8 nearest Grid Points:
		double iTemp = grid.i(xyzTemp[1]);
		double jTemp = grid.j(xyzTemp[2]);
		double kTemp = grid.k(xyzTemp[3]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;
		int k0 = std::floor(kTemp);	int k1 = k0 + 1;

		// Interpolate intensity from neighbouring 8 lattice points to temporary point:
		Inew[ijkd]
		= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 IntensityAt(grid.Index(i0,j0,k0), vTempIF), IntensityAt(grid.Index(i0,j0,k1), vTempIF), IntensityAt(grid.Index(i0,j1,k0), vTempIF), IntensityAt(grid.Index(i0,j1,k1), vTempIF),
		 IntensityAt(grid.Index(i1,j0,k0), vTempIF), IntensityAt(grid.Index(i1,j0,k1), vTempIF), IntensityAt(grid.Index(i1,j1,k0), vTempIF), IntensityAt(grid.Index(i1,j1,k1), vTempIF));
	}
	std::swap(I,Inew);
	std::swap(nx,nxNew);
	std::swap(ny,nyNew);
	std::swap(nz,nzNew);
}



void Radiation::StreamGeodesicStatic()
{
	exit_on_error("StreamGeodesicDynamic not supported yet.");
	//PROFILE_FUNCTION();
	//#pragma omp parallel for
	//for(int d=0; d<nDir; d++)
	//for(int j=2; j<grid.n2-2; j++)
	//for(int i=2; i<grid.n1-2; i++)
	//{
	//	int ijk = grid.Index(i,j,k);		// Index of lattice point ijk
	//	int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ijk
	//
	//	// Skip LPs which are inside BH:
	//	if(metric.InsideBH(i,j))
	//	{
	//		Inew[ijd] = 0;
	//		continue;
	//	}
	//
	//	// Geodesic Streaming:
	//	double s = 1;
	//	double alpha = metric.GetAlpha(ijk);
	//	Coordinate2<Coord> xTemp = grid.x12Coord(i,j);
    //    Tensor2<xy,IF> cxy = stencil.Cxy(d);
    //    Tensor2<Coord,IF> c = cxy.template Transform<Coord>(grid.xyCoord(i,j));
    //    Tensor3<Coord,IF> uIF(alpha, c[1]*alpha, c[2]*alpha);
    //    Tensor2<Coord,LF> vTempLF = Vec2ObservedByEulObs<Coord,IF,LF>(uIF,xTemp,metric);
	//	if(!metric.InsideBH(xTemp))
    //		s /= RK45_GeodesicEquation<-1>(grid.dt,xTemp,vTempLF,metric);
	//	Tensor2<Coord,IF> vTempIF = vTempLF.template Transform<IF>(metric.GetTetrad(xTemp));
	//	Tensor2<xy,IF> vTempIFxy = vTempIF.template Transform<xy>(xTemp);
	//
	//	// Skip temporary Grid Points inside BH:
	//	if(metric.InsideBH(xTemp))
	//	{
	//		Inew[ijd] = 0;
	//		continue;
	//	}
	//
	//	// Get 4 nearest Grid Points:
	//	double iTemp = grid.i(xTemp[1]);
	//	double jTemp = grid.j(xTemp[2]);
	//	int i0 = std::floor(iTemp);	int i1 = i0 + 1;
	//	int j0 = std::floor(jTemp);	int j1 = j0 + 1;
	//
	//	// Intensity interpolation:
	//	double angle = vTempIFxy.Angle();
	//	double intensityAt_i0j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0))) * IntensityAt(grid.Index(i0,j0),angle);
	//	double intensityAt_i0j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1))) * IntensityAt(grid.Index(i0,j1),angle);
	//	double intensityAt_i1j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0))) * IntensityAt(grid.Index(i1,j0),angle);
	//	double intensityAt_i1j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1))) * IntensityAt(grid.Index(i1,j1),angle);
	//
	//	// Interpolate intensity from neighbouring 4 lattice points to temporary point:
	//	Inew[ijd]
	//	= (s*s*s*s) * BilinearInterpolation(iTemp-i0,jTemp-j0,
	//	 intensityAt_i0j0,intensityAt_i0j1,
	//	 intensityAt_i1j0,intensityAt_i1j1);
	//}
	//// swap old with new array:
	//std::swap(I,Inew);
}
void Radiation::StreamGeodesicDynamic()
{
	exit_on_error("StreamGeodesicDynamic not supported yet.");
	//PROFILE_FUNCTION();
	//#pragma omp parallel for
	//for(int j=2; j<grid.n2-2; j++)
	//for(int i=2; i<grid.n1-2; i++)
	//{
	//	int ijk = grid.Index(i,j,k);
	//	Tensor2<Coord,IF> averageF = AverageF(i,j);
	//	constexpr double rotationThreshold = 1e-3;
	//	rotationNew[ijk] = (averageF.EuklNorm()>rotationThreshold) ? averageF.Angle() : rotation[ijk];
	//}
	//
	//#pragma omp parallel for
	//for(int d=0; d<nDir; d++)
	//for(int j=2; j<grid.n2-2; j++)
	//for(int i=2; i<grid.n1-2; i++)
	//{
	//	int ijk = grid.Index(i,j,k);		// Index of lattice point ijk
	//	int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ijk
	//
	//	// Skip LPs which are inside BH:
	//	if(metric.InsideBH(i,j))
	//	{
	//		Inew[ijd] = 0;
	//		continue;
	//	}
	//
	//	// Geodesic Streaming:
	//	double s = 1;
	//	double alpha = metric.GetAlpha(ijk);
	//	Coordinate2<Coord> xTemp = grid.x12Coord(i,j);
    //    Tensor2<xy,IF> cxy = stencil.Cxy(d,rotationNew[ijk]);
    //    Tensor2<Coord,IF> c = cxy.template Transform<Coord>(grid.xyCoord(i,j));
    //    Tensor3<Coord,IF> uIF(alpha, c[1]*alpha, c[2]*alpha);
    //    Tensor2<Coord,LF> vTempLF = Vec2ObservedByEulObs<Coord,IF,LF>(uIF,xTemp,metric);
	//	if(!metric.InsideBH(xTemp))
    //		s /= RK45_GeodesicEquation<-1>(grid.dt,xTemp,vTempLF,metric);
	//	Tensor2<Coord,IF> vTempIF = vTempLF.template Transform<IF>(metric.GetTetrad(xTemp));
	//	Tensor2<xy,IF> vTempIFxy = vTempIF.template Transform<xy>(xTemp);
	//
	//	// Skip temporary Grid Points inside BH:
	//	if(metric.InsideBH(xTemp))
	//	{
	//		Inew[ijd] = 0;
	//		continue;
	//	}
	//
	//	// Get 4 nearest Grid Points:
	//	double iTemp = grid.i(xTemp[1]);
	//	double jTemp = grid.j(xTemp[2]);
	//	int i0 = std::floor(iTemp);	int i1 = i0 + 1;
	//	int j0 = std::floor(jTemp);	int j1 = j0 + 1;
	//
	//	// Intensity interpolation:
	//	double angle = vTempIFxy.Angle();
	//	double intensityAt_i0j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0))) * IntensityAt(grid.Index(i0,j0),angle);
	//	double intensityAt_i0j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1))) * IntensityAt(grid.Index(i0,j1),angle);
	//	double intensityAt_i1j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0))) * IntensityAt(grid.Index(i1,j0),angle);
	//	double intensityAt_i1j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1))) * IntensityAt(grid.Index(i1,j1),angle);
	//
	//	// Interpolate intensity from neighbouring 4 lattice points to temporary point:
	//	Inew[ijd]
	//	= (s*s*s*s) * BilinearInterpolation(iTemp-i0,jTemp-j0,
	//	 intensityAt_i0j0,intensityAt_i0j1,
	//	 intensityAt_i1j0,intensityAt_i1j1);
	//}
	//// swap old with new arrays:
	//std::swap(I,Inew);
	//std::swap(rotation,rotationNew);
}



void Radiation::StreamCurvedStatic()
{
	exit_on_error("StreamCurvedStatic not supported yet.");
	//PROFILE_FUNCTION();
	//#pragma omp parallel for
	//for(int d0=0; d0<stencil.nTh; d0++)
	//for(int d1=0; d1<stencil.nPh; d1++)
	//for(int k=2; k<grid.nz-2; k++)
	//for(int j=2; j<grid.ny-2; j++)
	//for(int i=2; i<grid.nx-2; i++)
	//{
	//	int d = stencil.Index(d0,d1);	// Index of direction d
	//	int ijk = grid.Index(i,j,k);	// Index of lattice point ijk
	//	int ijkd = grid.Index(i,j,k,d);	// Index of population d at lattice point ijk
	//
	//	// Skip LPs which are inside BH:
	//	if(metric.InsideBH(grid.xyz(i,j,k)))
	//	{
	//		Inew[ijkd] = 0;
	//		continue;
	//	}
	//
	//	// Curved Fourier Streaming:
	//	double theta = stencil.Theta(d0, d1);
	//	double phi = stencil.Phi(d0, d1);
	//	double s = GetFrequencyShift(i, j, k, theta, phi);
	//	Coord xyzTemp = GetTempCoordinate(i, j, k, theta, phi);
	//
	//	// Skip temporary Grid Points inside BH:
	//	if(metric.InsideBH(xyzTemp))
	//	{
	//		Inew[ijkd] = 0;
	//		continue;
	//	}
	//
	//	Tensor4x4 inverseTetrad = metric.GetTetrad(xyzTemp).Invert();
	//	Tensor3 vTempLF = GetTemp3Velocity(i, j, k, theta, phi);
	//	Tensor3 vTempIF = TransformLFtoIF(vTempLF, inverseTetrad);
	//
	//	// Get 4 nearest Grid Points:
	//	double iTemp = grid.i(xyzTemp[1]);
	//	double jTemp = grid.j(xyzTemp[2]);
	//	double kTemp = grid.k(xyzTemp[3]);
	//	int i0 = std::floor(iTemp);	int i1 = i0 + 1;
	//	int j0 = std::floor(jTemp);	int j1 = j0 + 1;
	//	int k0 = std::floor(kTemp);	int k1 = k0 + 1;
	//
	//	// Intensity interpolation:
	//	double vTheta = vTempIF.Theta();
	//	double vPhi = vTempIF.Phi();
	//	double alpha = metric.GetAlpha(ijk);
	//	double intensityAt_i0j0k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k0))) * IntensityAt(grid.Index(i0,j0,k0),vTheta,vPhi);
	//	double intensityAt_i0j0k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k1))) * IntensityAt(grid.Index(i0,j0,k1),vTheta,vPhi);
	//	double intensityAt_i0j1k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k0))) * IntensityAt(grid.Index(i0,j1,k0),vTheta,vPhi);
	//	double intensityAt_i0j1k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k1))) * IntensityAt(grid.Index(i0,j1,k1),vTheta,vPhi);
	//	double intensityAt_i1j0k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k0))) * IntensityAt(grid.Index(i1,j0,k0),vTheta,vPhi);
	//	double intensityAt_i1j0k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k1))) * IntensityAt(grid.Index(i1,j0,k1),vTheta,vPhi);
	//	double intensityAt_i1j1k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k0))) * IntensityAt(grid.Index(i1,j1,k0),vTheta,vPhi);
	//	double intensityAt_i1j1k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k1))) * IntensityAt(grid.Index(i1,j1,k1),vTheta,vPhi);
	//
	//	// Interpolate intensity from neighbouring 4 lattice points to temporary point:
	//	Inew[ijkd]
	//	= (s*s*s*s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
	//	 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
	//	 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
	//}
	//// swap old with new array:
	//std::swap(I,Inew);
}
void Radiation::StreamCurvedDynamic()
{
	exit_on_error("StreamCurvedDynamic not supported yet.");
	// TODO: write this code

	//PROFILE_FUNCTION();
	//#pragma omp parallel for
	//for(int j=2; j<grid.n2-2; j++)
	//for(int i=2; i<grid.n1-2; i++)
	//{
	//	int ijk = grid.Index(i,j,k);
	//	Tensor2<Coord,IF> averageF = AverageF(i,j);
	//	constexpr double rotationThreshold = 1e-3;
	//	rotationNew[ijk] = (averageF.EuklNorm()>rotationThreshold) ? averageF.Angle() : rotation[ijk];
	//}
	//
	//#pragma omp parallel for
	//for(int d=0; d<nDir; d++)
	//for(int j=2; j<grid.n2-2; j++)
	//for(int i=2; i<grid.n1-2; i++)
	//{
	//	int ijk = grid.Index(i,j,k);		// Index of lattice point ijk
	//	int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ijk
	//
	//	// Skip LPs which are inside BH:
	//	if(metric.InsideBH(i,j))
	//	{
	//		Inew[ijd] = 0;
	//		continue;
	//	}
	//
	//	// Curved Fourier Streaming:
	//	double phi = stencil.Phi(d,rotationNew[ijk]);
	//	double s = GetFrequencyShift(i,j,phi);
	//	Coordinate2<Coord> xTemp = GetTempCoordinate(i,j,phi);
	//	Tensor2<Coord,LF> vTempLF = GetTemp2Velocity(i,j,phi);
	//	Tensor2<Coord,IF> vTempIF = vTempLF.template Transform<IF>(metric.GetTetrad(xTemp));
	//	Tensor2<xy,IF> vTempIFxy = vTempIF.template Transform<xy>(xTemp);
	//
	//	// Skip temporary Grid Points inside BH:
	//	if(metric.InsideBH(xTemp))
	//	{
	//		Inew[ijd] = 0;
	//		continue;
	//	}
	//
	//	// Get 4 nearest Grid Points:
	//	double iTemp = grid.i(xTemp[1]);
	//	double jTemp = grid.j(xTemp[2]);
	//	int i0 = std::floor(iTemp);	int i1 = i0 + 1;
	//	int j0 = std::floor(jTemp);	int j1 = j0 + 1;
	//	int i0j0 = grid.Index(i0,j0);
	//	int i0j1 = grid.Index(i0,j1);
	//	int i1j0 = grid.Index(i1,j0);
	//	int i1j1 = grid.Index(i1,j1);
	//
	//	// TEST:
	//	// double angle00 = ((vTempLF.template Transform<IF>(metric.GetTetrad(i0j0))).template Transform<xy>(grid.x12Coord(i0,j0))).Angle();
	//	// double angle01 = ((vTempLF.template Transform<IF>(metric.GetTetrad(i0j1))).template Transform<xy>(grid.x12Coord(i0,j1))).Angle();
	//	// double angle10 = ((vTempLF.template Transform<IF>(metric.GetTetrad(i1j0))).template Transform<xy>(grid.x12Coord(i1,j0))).Angle();
	//	// double angle11 = ((vTempLF.template Transform<IF>(metric.GetTetrad(i1j1))).template Transform<xy>(grid.x12Coord(i1,j1))).Angle();
	//
	//	// Intensity interpolation:
	//	double angle = vTempIFxy.Angle();
	//	double alpha = metric.GetAlpha(ijk);
	//	double intensityAt_i0j0 = MyPow<4>(alpha / metric.GetAlpha(i0j0)) * IntensityAt(i0j0,angle);
	//	double intensityAt_i0j1 = MyPow<4>(alpha / metric.GetAlpha(i0j1)) * IntensityAt(i0j1,angle);
	//	double intensityAt_i1j0 = MyPow<4>(alpha / metric.GetAlpha(i1j0)) * IntensityAt(i1j0,angle);
	//	double intensityAt_i1j1 = MyPow<4>(alpha / metric.GetAlpha(i1j1)) * IntensityAt(i1j1,angle);
	//
	//	// Interpolate intensity from neighbouring 4 lattice points to temporary point:
	//	Inew[ijd]
	//	= (s*s*s*s) * BilinearInterpolation(iTemp-i0,jTemp-j0,
	//	 intensityAt_i0j0,intensityAt_i0j1,
	//	 intensityAt_i1j0,intensityAt_i1j1);
	//}
	//// swap old with new arrays:
	//std::swap(I,Inew);
	//std::swap(rotation,rotationNew);
}



void Radiation::Collide()
{
	PROFILE_FUNCTION();
	// TODO: Steife DGL?
	
	#pragma omp parallel for
	for(int k = 0; k < metric.grid.nz; k++)
	for(int j = 0; j < metric.grid.ny; j++)
	for(int i = 0; i < metric.grid.nx; i++)
	{
		if(metric.InsideBH(grid.xyz(i,j,k)))
			continue;

		int ijk = grid.Index(i,j,k);

		// Simulate stationary fluid, u^k=(0,0):
		double alpha = metric.GetAlpha(ijk);
		Tensor3 u(0.0);	// fluid 3 velocity as seen by Eulerian observer
		double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ijk)));	// Lorentz factor
		double uDotF = Fx_LF[ijk] * u[1] + Fy_LF[ijk] * u[2] + Fz_LF[ijk] * u[3];	// F^i u_i
		double uuDotP =
		Pxx_LF[ijk] * u[1] * u[1] + Pyy_LF[ijk] * u[2] * u[2] + Pzz_LF[ijk] * u[3] * u[3]
		+ 2.0 * (Pxy_LF[ijk] * u[1] * u[2] + Pxz_LF[ijk] * u[1] * u[3] + Pyz_LF[ijk] * u[2] * u[3]) ;	// P^ijk u_i u_j
		double fluidE = W * W * (E_LF[ijk]  - 2.0 * uDotF + uuDotP);
		
		for(int d1 = 0; d1 < stencil.nPh; d1++)
		for(int d0 = 0; d0 < stencil.nTh; d0++)
		{
			int d = stencil.Index(d0,d1);
			int ijd = ijk + d*metric.grid.nxyz;
			double A = W * (1.0 - Tensor3::Dot(stencil.Cxyz(d0,d1), u));

			double Gamma = stencil.W(d0,d1) * (eta[ijk] + kappa0[ijk]*fluidE) / (A*A*A) - A*I[ijd] * (kappaA[ijk] + kappa0[ijk]);
			I[ijd] += alpha * metric.grid.dt * Gamma;
		}
	}
}



void Radiation::RunSimulation(Config config)
{
	// -------------------- Initialization --------------------
	int timeSteps = ceil(config.simTime / grid.dt);
	config.simTime = timeSteps * grid.dt;
	Log logger(config.name, config.simTime, stencil, lebedevStencil, metric);

	NormalizeInitialDirections();
	LoadInitialData();
	ComputeMomentsIF();
	ComputeMomentsLF();
	NormalizeInitialData();
	LoadInitialData(); // loads normalized data
	UpdateSphericalHarmonicsCoefficients();

	// Initial data output:
	if (config.printToTerminal)
	{
		std::cout << " sigma        = " << stencil.sigma << "\n";
		std::cout << " nx           = " << grid.nx << "\n";
		std::cout << " ny           = " << grid.ny << "\n";
		std::cout << " nz           = " << grid.nz << "\n";
		std::cout << " nTh          = " << stencil.nTh << "\n";
		std::cout << " nPh          = " << stencil.nPh << "\n";
		std::cout << " nSphHarm     = " << lebedevStencil.nDir << "\n";
		std::cout << " simTime      = " << config.simTime << "\n";
		std::cout << " dt           = " << grid.dt << "\n";
		std::cout << " timeSteps    = " << timeSteps << "\n";
		std::cout << " filesToWrite = " << timeSteps / config.writeFrequency << std::endl;
	}
	// --------------------------------------------------------



	Profiler::Session& session = Profiler::Session::Get();
	session.Start(config.name, "output/" + config.name + "/profileResults.json");
	// ----------------- Main simulation Loop -----------------
	{
		PROFILE_SCOPE("Total Time");
		for(int n=0; n<timeSteps; n++)
		{
			if (config.printToTerminal)
			{ std::cout << n << "," << std::flush; }
			ComputeMomentsIF();
			Collide();

			if (config.writeData && (n % config.writeFrequency) == 0)
			{
				ComputeMomentsLF();
				grid.WriteFrametoCsv(n*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, n, logger.directoryPath + "/Moments", config.name);
			}

			if (config.updateSphericalHarmonics)
			{ UpdateSphericalHarmonicsCoefficients(); }

			SetIntensitiesNorthSouth();

			// Streaming:
			switch (streamingType)
			{
				case(StreamingType::FlatStatic):		StreamFlatStatic();		 break;
				case(StreamingType::FlatDynamic):		StreamFlatDynamic();	 break;
				case(StreamingType::GeodesicStatic):	StreamGeodesicStatic();	 break;
				case(StreamingType::GeodesicDynamic):	StreamGeodesicDynamic(); break;
				case(StreamingType::CurvedStatic):		StreamCurvedStatic();	 break;
				case(StreamingType::CurvedDynamic):		StreamCurvedDynamic();	 break;
			}

			if (config.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if (config.writeData)
		{
			ComputeMomentsLF();
			grid.WriteFrametoCsv(timeSteps*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, timeSteps, logger.directoryPath + "/Moments", config.name);
		}
	}
	// --------------------------------------------------------
	session.End();



	// ---------------------- Termination ---------------------
	std::vector<std::string> names = session.GetAllFunctionNames();
	if (config.printToTerminal)
	{ std::cout << std::endl; }
	for(int i=0; i<names.size(); i++)
	{
		if (config.printToTerminal)
			session.PrintFunctionDuration(names[i]);
		logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
	}
	// --------------------------------------------------------
}