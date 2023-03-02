#include "Radiation.h"

std::string StreamingName(int n)
{
	std::string name("unknown");
	switch (n)
	{
   		case 0: { name = "FlatStatic";		} break;
   		case 1: { name = "FlatDynamic";		} break;
   		case 2: { name = "CurvedStatic";	} break;
   		case 3: { name = "CurvedDynamic";	} break;
		default: { exit_on_error("Invalid StreamingType"); }
	}
	return name;
}



Radiation::Radiation(Metric& metric, Stencil& stencil, LebedevStencil& lebedevStencil, Camera& camera, StreamingType streamingType):
grid(metric.grid), metric(metric), stencil(stencil), lebedevStencil(lebedevStencil), camera(camera), streamingType(streamingType)
{
	isInitialGridPoint = new bool[grid.nxyz]();
	initialE.resize(grid.nxyz);
	initialNx.resize(grid.nxyz);
	initialNy.resize(grid.nxyz);
	initialNz.resize(grid.nxyz);
	initialKappa0.resize(grid.nxyz);
	initialKappa1.resize(grid.nxyz);
	initialKappaA.resize(grid.nxyz);
	initialEta.resize(grid.nxyz);
	nx.resize(grid.nxyz);
	ny.resize(grid.nxyz);
	nz.resize(grid.nxyz);
	nxNew.resize(grid.nxyz);
	nyNew.resize(grid.nxyz);
	nzNew.resize(grid.nxyz);
	E.resize(grid.nxyz);
	Fx.resize(grid.nxyz);
	Fy.resize(grid.nxyz);
	Fz.resize(grid.nxyz);
	Pxx.resize(grid.nxyz);
	Pxy.resize(grid.nxyz);
	Pxz.resize(grid.nxyz);
	Pyy.resize(grid.nxyz);
	Pyz.resize(grid.nxyz);
	Pzz.resize(grid.nxyz);
	E_LF .resize(grid.nxyz);
	Fx_LF.resize(grid.nxyz);
	Fy_LF.resize(grid.nxyz);
	Fz_LF.resize(grid.nxyz);
	Pxx_LF.resize(grid.nxyz);
	Pxy_LF.resize(grid.nxyz);
	Pxz_LF.resize(grid.nxyz);
	Pyy_LF.resize(grid.nxyz);
	Pyz_LF.resize(grid.nxyz);
	Pzz_LF.resize(grid.nxyz);
	kappa0.resize(grid.nxyz);
	kappa1.resize(grid.nxyz);
	kappaA.resize(grid.nxyz);
	eta.resize(grid.nxyz);
	I.resize((size_t)grid.nxyz * (size_t)stencil.nThPh);
	Inew.resize((size_t)grid.nxyz * (size_t)stencil.nThPh);
	Inorth.resize(grid.nxyz);
	Isouth.resize(grid.nxyz);
	coefficientsS.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsX.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsY.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsZ.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsCx.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsCy.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsCz.resize(grid.nxyz * lebedevStencil.nCoefficients);

	// Initialize all stencil directions to north pole:
	PARALLEL_FOR(1)
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		nx[ijk] = 0;
		ny[ijk] = 0;
		nz[ijk] = 1;
	}
}



Radiation::~Radiation()
{
	delete[] isInitialGridPoint;
}




int Radiation::Index(int ijk, int d)
{
	return ijk + d * grid.nxyz;
}
int Radiation::Index(int ijk, int d0, int d1)
{
	int d = stencil.Index(d0,d1);
	return ijk + d * grid.nxyz;
}
int Radiation::Index(int i, int j, int k, int d)
{
	int ijk = grid.Index(i,j,k);
	return ijk + d * grid.nxyz;
}
int Radiation::Index(int i, int j, int k, int d0, int d1)
{
	int ijk = grid.Index(i,j,k);
	int d = stencil.Index(d0,d1);
	return ijk + d * grid.nxyz;
}

int Radiation::HarmonicIndex(int f, int ijk)
{
	// faster:
	return f + ijk * lebedevStencil.nCoefficients;
	// slower:
	// return f * grid.nxyz + ijk;
}



void Radiation::NormalizeInitialDirections()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(1)
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
	DEBUG_FUNCTION();
	bool isDynamicStreaming = (streamingType == StreamingType::FlatDynamic || streamingType == StreamingType::CurvedDynamic);
	srand((unsigned) time(NULL));

	PARALLEL_FOR(1)
	for(int ijk=0; ijk<grid.nxyz; ijk++)
	{
		kappa0[ijk] = initialKappa0[ijk];
		kappa1[ijk] = initialKappa1[ijk];
		kappaA[ijk] = initialKappaA[ijk];
		eta[ijk] = initialEta[ijk];
		
		if(isInitialGridPoint[ijk])
		{
			nx[ijk] = initialNx[ijk];
			ny[ijk] = initialNy[ijk];
			nz[ijk] = initialNz[ijk];

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
					I[Index(ijk,d)] = initialE[ijk];
				Inorth[ijk] = initialE[ijk];
				Isouth[ijk] = initialE[ijk];
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
					I[Index(ijk,d)] = initialE[ijk] * exp(stencil.sigma * Tensor3::Dot(n, p));
				}
				Inorth[ijk] = initialE[ijk] * exp(stencil.sigma * Tensor3::Dot(n, Tensor3(0,0, 1)));
				Isouth[ijk] = initialE[ijk] * exp(stencil.sigma * Tensor3::Dot(n, Tensor3(0,0,-1)));
			}
		}
	}
}



void Radiation::NormalizeInitialIntensities()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(1)
	for(int ijk=0; ijk<grid.nxyz; ijk++)
		if(isInitialGridPoint[ijk])
			initialE[ijk] *= initialE[ijk] / E_LF[ijk];
}



void Radiation::UpdateSphericalHarmonicsCoefficients()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(3)
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
    	    	s *= RK45_GeodesicEquation<-1>(grid.dt, xyz, v, metric);
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
			int fijk = HarmonicIndex(f,ijk);
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
	DEBUG_FUNCTION();
	constexpr double fourPiInv = 1.0 / (4.0 * M_PI);
	
	PARALLEL_FOR(1)
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
			int index = Index(ijk,d);
			double c = stencil.W(d0,d1) * I[index];
			E[ijk]   += c;
			Fx[ijk]  += c * cxyz[1];
			Fy[ijk]  += c * cxyz[2];
			Fz[ijk]  += c * cxyz[3];
			Pxx[ijk] += c * cxyz[1] * cxyz[1];
			Pxy[ijk] += c * cxyz[1] * cxyz[2];
			Pxz[ijk] += c * cxyz[1] * cxyz[3];
			Pyy[ijk] += c * cxyz[2] * cxyz[2];
			Pyz[ijk] += c * cxyz[2] * cxyz[3];
			Pzz[ijk] += c * cxyz[3] * cxyz[3];
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
	DEBUG_FUNCTION();
	PARALLEL_FOR(1)
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



Coord Radiation::GetTempCoordinate(int ijk, Tensor3 direction)
{
	Coord xyzTemp;
	double evaluationCoefficients[lebedevStencil.nCoefficients];

	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsX[HarmonicIndex(f,ijk)];
	xyzTemp[1] = SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsY[HarmonicIndex(f,ijk)];
	xyzTemp[2] = SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsZ[HarmonicIndex(f,ijk)];
	xyzTemp[3] = SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);

	return xyzTemp;
}
Tensor3 Radiation::GetTemp3Velocity(int ijk, Tensor3 direction)
{
	Tensor3 vTempIF;
	double evaluationCoefficients[lebedevStencil.nCoefficients];

	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsCx[HarmonicIndex(f,ijk)];
	vTempIF[1] = SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsCy[HarmonicIndex(f,ijk)];
	vTempIF[2] = SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsCz[HarmonicIndex(f,ijk)];
	vTempIF[3] = SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);

	return vTempIF;
}
double Radiation::GetFrequencyShift(int ijk, Tensor3 direction)
{
	double evaluationCoefficients[lebedevStencil.nCoefficients];
	for(int f=0; f<lebedevStencil.nCoefficients; f++)
		evaluationCoefficients[f] = coefficientsS[HarmonicIndex(f,ijk)];
	return SphericalHarmonics::GetValue(direction,evaluationCoefficients,lebedevStencil.nCoefficients);
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
    double I00 = (th0 == -1)          ? Inorth[ijk] : I[Index(ijk,th0,ph0)];
    double I01 = (th0 == -1)          ? Inorth[ijk] : I[Index(ijk,th0,ph1)];
    double I10 = (th1 == stencil.nTh) ? Isouth[ijk] : I[Index(ijk,th1,ph0)];
    double I11 = (th1 == stencil.nTh) ? Isouth[ijk] : I[Index(ijk,th1,ph1)];

	// Fractional part of sphere grid index. Near poles it must be streched:
	// North: from [0.5, 1.0] to [0.0, 1.0]
	// South: from [0.0, 0.5] to [0.0, 1.0]
	double tfrac = th - floor(th);
	double pfrac = ph - floor(ph);
    if (th0 == -1)					// North
        tfrac = 2.0 * (tfrac - 0.5);
    else if (th1 == stencil.nTh)	// South
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



void Radiation::GetNewRotation()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
	{
		int ijk = grid.Index(i,j,k);
		Tensor3 averageF = AverageF(i,j,k);
		Tensor3 n = (averageF.EuklNorm()>normThreshhold) ? averageF.EuklNormalized() : Tensor3(nx[ijk],ny[ijk],nz[ijk]);
		nxNew[ijk] = n[1];
		nyNew[ijk] = n[2];
		nzNew[ijk] = n[3];
	}
}



void Radiation::StreamFlatStatic()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(5)
	for(int d1=0; d1<stencil.nPh; d1++)
	for(int d0=0; d0<stencil.nTh; d0++)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamFlatKernal<Bulk,Static>(i,j,k,d0,d1);

	std::swap(I,Inew);
}
void Radiation::StreamFlatDynamic()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(5)
	for(int d1=0; d1<stencil.nPh; d1++)
	for(int d0=0; d0<stencil.nTh; d0++)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamFlatKernal<Bulk,Dynamic>(i,j,k,d0,d1);

	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamFlatKernal<North,Dynamic>(i,j,k,0,0);
			
	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamFlatKernal<South,Dynamic>(i,j,k,0,0);

	std::swap(I,Inew);
	std::swap(nx,nxNew);
	std::swap(ny,nyNew);
	std::swap(nz,nzNew);
}
template<class IntensityType, class StaticOrDynamic>
void Radiation::StreamFlatKernal(int i, int j, int k, int d0, int d1)
{
	int ijk = grid.Index(i,j,k);	// Index of lattice point ijk
	int d = stencil.Index(d0,d1);	// Index of direction d
	int index = Index(ijk,d);		// Index of population d at lattice point ijk

	// Flat Streaming:
	glm::vec3 from(0,0,1);
	glm::vec3 to(nxNew[ijk],nyNew[ijk],nzNew[ijk]);
	glm::quat q(from,to);

	// Get temp velocity:
	Tensor3 vTempIF;
	if constexpr(std::is_same<StaticOrDynamic,Dynamic>::value)
	{
		if constexpr(std::is_same<IntensityType,Bulk>::value)
			vTempIF = q * stencil.Cxyz(d0,d1);
		if constexpr(std::is_same<IntensityType,North>::value)
			vTempIF = q * Tensor3(0,0,1);
		if constexpr(std::is_same<IntensityType,South>::value)
			vTempIF = q * Tensor3(0,0,-1);
	}
	if constexpr(std::is_same<StaticOrDynamic,Static>::value)
		vTempIF = stencil.Cxyz(d0,d1);

	// Get temp lattice point:
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
	if constexpr(std::is_same<StaticOrDynamic,Dynamic>::value)
	{
		if constexpr(std::is_same<IntensityType,Bulk>::value)
			Inew[index]
			= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
			 IntensityAt(grid.Index(i0,j0,k0), vTempIF), IntensityAt(grid.Index(i0,j0,k1), vTempIF), IntensityAt(grid.Index(i0,j1,k0), vTempIF), IntensityAt(grid.Index(i0,j1,k1), vTempIF),
			 IntensityAt(grid.Index(i1,j0,k0), vTempIF), IntensityAt(grid.Index(i1,j0,k1), vTempIF), IntensityAt(grid.Index(i1,j1,k0), vTempIF), IntensityAt(grid.Index(i1,j1,k1), vTempIF));
		if constexpr(std::is_same<IntensityType,North>::value)
			Inorth[ijk]
			= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
			 IntensityAt(grid.Index(i0,j0,k0), vTempIF), IntensityAt(grid.Index(i0,j0,k1), vTempIF), IntensityAt(grid.Index(i0,j1,k0), vTempIF), IntensityAt(grid.Index(i0,j1,k1), vTempIF),
			 IntensityAt(grid.Index(i1,j0,k0), vTempIF), IntensityAt(grid.Index(i1,j0,k1), vTempIF), IntensityAt(grid.Index(i1,j1,k0), vTempIF), IntensityAt(grid.Index(i1,j1,k1), vTempIF));
		if constexpr(std::is_same<IntensityType,South>::value)
			Isouth[ijk]
			= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
			 IntensityAt(grid.Index(i0,j0,k0), vTempIF), IntensityAt(grid.Index(i0,j0,k1), vTempIF), IntensityAt(grid.Index(i0,j1,k0), vTempIF), IntensityAt(grid.Index(i0,j1,k1), vTempIF),
			 IntensityAt(grid.Index(i1,j0,k0), vTempIF), IntensityAt(grid.Index(i1,j0,k1), vTempIF), IntensityAt(grid.Index(i1,j1,k0), vTempIF), IntensityAt(grid.Index(i1,j1,k1), vTempIF));	
	}
	if constexpr(std::is_same<StaticOrDynamic,Static>::value)
		Inew[index]
		= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		I[Index(i0,j0,k0,d)], I[Index(i0,j0,k1,d)], I[Index(i0,j1,k0,d)], I[Index(i0,j1,k1,d)],
		I[Index(i1,j0,k0,d)], I[Index(i1,j0,k1,d)], I[Index(i1,j1,k0,d)], I[Index(i1,j1,k1,d)]);
}



void Radiation::StreamCurvedStatic()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(5)
	for(int d1=0; d1<stencil.nPh; d1++)
	for(int d0=0; d0<stencil.nTh; d0++)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamCurvedKernal<Bulk,Static>(i,j,k,d0,d1);

	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamCurvedKernal<North,Dynamic>(i,j,k,0,0);
			
	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamCurvedKernal<South,Dynamic>(i,j,k,0,0);

	std::swap(I,Inew);
}
void Radiation::StreamCurvedDynamic()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	PARALLEL_FOR(5)
	for(int d1=0; d1<stencil.nPh; d1++)
	for(int d0=0; d0<stencil.nTh; d0++)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamCurvedKernal<Bulk,Dynamic>(i,j,k,d0,d1);
		
	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamCurvedKernal<North,Dynamic>(i,j,k,0,0);
			
	PARALLEL_FOR(3)
	for(int k=2; k<grid.nz-2; k++)
	for(int j=2; j<grid.ny-2; j++)
	for(int i=2; i<grid.nx-2; i++)
		StreamCurvedKernal<South,Dynamic>(i,j,k,0,0);

	std::swap(I,Inew);
	std::swap(nx,nxNew);
	std::swap(ny,nyNew);
	std::swap(nz,nzNew);
}
template<class IntensityType, class StaticOrDynamic>
void Radiation::StreamCurvedKernal(int i, int j, int k, int d0, int d1)
{
	int ijk = grid.Index(i,j,k);	// Index of lattice point ijk
	int d = stencil.Index(d0,d1);	// Index of direction d
	int index = Index(ijk,d);		// Index of population d at lattice point ijk

	// Skip LPs which are inside BH:
	if(metric.InsideBH(grid.xyz(i,j,k)))
	{
		Inew[index] = 0;
		return;
	}

	// Get velocity direction in IF:
	Tensor3 direction;
	if constexpr(std::is_same<StaticOrDynamic,Dynamic>::value)
	{
		// Prepare stencil rotation via quaternion:
		glm::vec3 from(0,0,1);
		glm::vec3 to(nxNew[ijk],nyNew[ijk],nzNew[ijk]);
		glm::quat q(from,to);

		// Get temp velocity in local inertial frame:
		if constexpr(std::is_same<IntensityType,Bulk>::value)
			direction = q * stencil.Cxyz(d0,d1);
		if constexpr(std::is_same<IntensityType,North>::value)
			direction = q * Tensor3(0,0,1);
		if constexpr(std::is_same<IntensityType,South>::value)
			direction = q * Tensor3(0,0,-1);
	}
	if constexpr(std::is_same<StaticOrDynamic,Static>::value)
	{
		if constexpr(std::is_same<IntensityType,Bulk>::value)
			direction = stencil.Cxyz(d0,d1);
		if constexpr(std::is_same<IntensityType,North>::value)
			direction = Tensor3(0,0,1);
		if constexpr(std::is_same<IntensityType,South>::value)
			direction = Tensor3(0,0,-1);
	}

	// Frequency Shift (s), Coordinates (xyz), and Velocity (v) at temp lattice point:
	double s = GetFrequencyShift(ijk, direction);
	Coord xyzTemp = GetTempCoordinate(ijk, direction);
	Tensor4x4 inverseTetrad = metric.GetTetrad(xyzTemp).Invert();
	Tensor3 vTempLF = GetTemp3Velocity(ijk, direction);
	Tensor3 vTempIF = TransformLFtoIF(vTempLF, inverseTetrad);

	// Skip temporary Grid Points inside BH:
	if(metric.InsideBH(xyzTemp))
	{
		Inew[index] = 0;
		return;
	}

	// Get 8 nearest Grid Points:
	double iTemp = grid.i(xyzTemp[1]);
	double jTemp = grid.j(xyzTemp[2]);
	double kTemp = grid.k(xyzTemp[3]);
	int i0 = std::floor(iTemp);	int i1 = i0 + 1;
	int j0 = std::floor(jTemp);	int j1 = j0 + 1;
	int k0 = std::floor(kTemp);	int k1 = k0 + 1;

	// Intensity interpolation:
	double alpha = metric.GetAlpha(ijk);
	double intensityAt_i0j0k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k0))) * IntensityAt(grid.Index(i0,j0,k0),vTempIF);
	double intensityAt_i0j0k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k1))) * IntensityAt(grid.Index(i0,j0,k1),vTempIF);
	double intensityAt_i0j1k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k0))) * IntensityAt(grid.Index(i0,j1,k0),vTempIF);
	double intensityAt_i0j1k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k1))) * IntensityAt(grid.Index(i0,j1,k1),vTempIF);
	double intensityAt_i1j0k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k0))) * IntensityAt(grid.Index(i1,j0,k0),vTempIF);
	double intensityAt_i1j0k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k1))) * IntensityAt(grid.Index(i1,j0,k1),vTempIF);
	double intensityAt_i1j1k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k0))) * IntensityAt(grid.Index(i1,j1,k0),vTempIF);
	double intensityAt_i1j1k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k1))) * IntensityAt(grid.Index(i1,j1,k1),vTempIF);

	// Interpolate intensity from neighbouring 4 lattice points to temporary point:
	if constexpr(std::is_same<IntensityType,Bulk>::value)
		Inew[index]
		= MyPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
		 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
	if constexpr(std::is_same<IntensityType,North>::value)
		Inorth[index]
		= MyPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
		 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
	if constexpr(std::is_same<IntensityType,South>::value)
		Isouth[index]
		= MyPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
		 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
}



void Radiation::Collide()
{
	PROFILE_FUNCTION();
	DEBUG_FUNCTION();
	// TODO: Steife DGL?
	
	PARALLEL_FOR(3)
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
		+ 2.0 * (Pxy_LF[ijk] * u[1] * u[2] + Pxz_LF[ijk] * u[1] * u[3] + Pyz_LF[ijk] * u[2] * u[3]);	// P^ij u_i u_j
		double fluidE = W * W * (E_LF[ijk] - 2.0 * uDotF + uuDotP);
		
		for(int d1 = 0; d1 < stencil.nPh; d1++)
		for(int d0 = 0; d0 < stencil.nTh; d0++)
		{
			int d = stencil.Index(d0,d1);
			int index = Index(ijk,d);
			double A = W * (1.0 - Tensor3::Dot(stencil.Cxyz(d0,d1), u));

			double Gamma = stencil.W(d0,d1) * (eta[ijk] + kappa0[ijk]*fluidE) / (A*A*A) - A*I[index] * (kappaA[ijk] + kappa0[ijk]);
			I[index] += alpha * metric.grid.dt * Gamma;
		}
	}
}



void Radiation::TakePicture()
{
	#pragma omp parallel for
	for(int ij=0; ij<camera.pixelCount; ij++)
	{
		Coord pixel = camera.xyz(ij);
		if(grid.OutsideDomain(pixel))
		{
			camera.image[ij] = 0;
			continue;
		}

		Tensor4x4 inverseTetrad = metric.GetTetrad(pixel).Invert();
		Tensor3 vTempLF = camera.orthogonalPassThrough;
		Tensor3 vTempIF = TransformLFtoIF(vTempLF, inverseTetrad);

		// Get 8 nearest Grid Points:
		double iTemp = grid.i(pixel[1]);
		double jTemp = grid.j(pixel[2]);
		double kTemp = grid.k(pixel[3]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;
		int k0 = std::floor(kTemp);	int k1 = k0 + 1;

		// Intensity interpolation:
		double alpha = metric.GetAlpha(pixel);
		double intensityAt_i0j0k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k0))) * IntensityAt(grid.Index(i0,j0,k0),vTempIF);
		double intensityAt_i0j0k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k1))) * IntensityAt(grid.Index(i0,j0,k1),vTempIF);
		double intensityAt_i0j1k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k0))) * IntensityAt(grid.Index(i0,j1,k0),vTempIF);
		double intensityAt_i0j1k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k1))) * IntensityAt(grid.Index(i0,j1,k1),vTempIF);
		double intensityAt_i1j0k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k0))) * IntensityAt(grid.Index(i1,j0,k0),vTempIF);
		double intensityAt_i1j0k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k1))) * IntensityAt(grid.Index(i1,j0,k1),vTempIF);
		double intensityAt_i1j1k0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k0))) * IntensityAt(grid.Index(i1,j1,k0),vTempIF);
		double intensityAt_i1j1k1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k1))) * IntensityAt(grid.Index(i1,j1,k1),vTempIF);	
		
		camera.image[ij]
		= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
		 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
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
	NormalizeInitialIntensities();
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
		std::cout << " dx           = " << grid.dx << "\n";
		std::cout << " dy           = " << grid.dy << "\n";
		std::cout << " dz           = " << grid.dz << "\n";
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
			if(config.printToTerminal)
			{ std::cout << n << "," << std::flush; }
			ComputeMomentsIF();
			Collide();

			if(config.writeData && (n % config.writeFrequency) == 0)
			{
				ComputeMomentsLF();
				grid.WriteFrametoCsv(n*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, n, logger.directoryPath + "/Moments", config.name);
			}
			if(config.useCamera && (n % config.writeFrequency) == 0)
			{
				TakePicture();
				camera.WriteImagetoCsv(n*grid.dt, n, logger.directoryPath + "/CameraImages");
			}

			if (config.updateSphericalHarmonics)
			{ UpdateSphericalHarmonicsCoefficients(); }

			// Streaming:
			switch(streamingType)
			{
				case(StreamingType::FlatStatic):		StreamFlatStatic();							break;
				case(StreamingType::FlatDynamic):		GetNewRotation(); StreamFlatDynamic();		break;
				case(StreamingType::CurvedStatic):		StreamCurvedStatic();	 					break;
				case(StreamingType::CurvedDynamic):		GetNewRotation(); StreamCurvedDynamic();	break;
			}

			if(config.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if(config.writeData)
		{
			ComputeMomentsLF();
			grid.WriteFrametoCsv(timeSteps*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, timeSteps, logger.directoryPath + "/Moments", config.name);
		}
		if(config.useCamera)
		{
			TakePicture();
			camera.WriteImagetoCsv(timeSteps*grid.dt, timeSteps, logger.directoryPath + "/CameraImages");
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