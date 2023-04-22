#include "Radiation.h"

/*

Radiation::Radiation(Metric& metric, MyStencil& stencil, LebedevStencil& lebedevStencil, Camera& camera, StreamingType streamingType):
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
	q.resize(grid.nxyz);
	qNew.resize(grid.nxyz);
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
	I.resize(grid.nxyz * stencil.nDir);
	Inew.resize(grid.nxyz * stencil.nDir);
	Inorth.resize(grid.nxyz * stencil.nDir);
	Isouth.resize(grid.nxyz * stencil.nDir);
	coefficientsS.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsX.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsY.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsZ.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsCx.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsCy.resize(grid.nxyz * lebedevStencil.nCoefficients);
	coefficientsCz.resize(grid.nxyz * lebedevStencil.nCoefficients);

	// Initialize all Quaternions to identity:
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
	{
		q[ijk] = glm::quat(1,0,0,0);
		qNew[ijk] = glm::quat(1,0,0,0);
	}
}
Radiation::~Radiation()
{
	delete[] isInitialGridPoint;
}




size_t Radiation::Index(size_t ijk, size_t d)
{
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return ijk * stencil.nDir + d;
	#endif
}
size_t Radiation::Index(size_t ijk, size_t d0, size_t d1)
{
	size_t d = stencil.Index(d0,d1);
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return ijk * stencil.nDir + d;
	#endif
}
size_t Radiation::Index(size_t i, size_t j, size_t k, size_t d)
{
	size_t ijk = grid.Index(i,j,k);
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return ijk * stencil.nDir + d;
	#endif
}
size_t Radiation::Index(size_t i, size_t j, size_t k, size_t d0, size_t d1)
{
	size_t ijk = grid.Index(i,j,k);
	size_t d = stencil.Index(d0,d1);
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return ijk * stencil.nDir + d;
	#endif
}
size_t Radiation::HarmonicIndex(size_t f, size_t ijk)
{
	return f + ijk * lebedevStencil.nCoefficients;
}



void Radiation::NormalizeInitialDirections()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
	{
		Tensor3 n(initialNx[ijk], initialNy[ijk], initialNz[ijk]);
		double norm = n.EuklNorm();

		// If norm to low we assume no direction.
		if(norm < 0.0001)
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
	bool isDynamicStreaming = (streamingType == StreamingType::FlatDynamic || streamingType == StreamingType::CurvedDynamic);
	srand((unsigned) time(NULL));

	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
	{
		kappa0[ijk] = initialKappa0[ijk];
		kappa1[ijk] = initialKappa1[ijk];
		kappaA[ijk] = initialKappaA[ijk];
		eta[ijk] = initialEta[ijk];
		
		if(isInitialGridPoint[ijk])
		{
			if(initialNx[ijk] == 0 && initialNy[ijk] == 0 && initialNz[ijk] == 0)
			{// Uniform intensity distribution:
				if (isDynamicStreaming)
				{// Random direction for uniformity:
					Tensor3 n(0.0);
					n[1] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
					n[2] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
					n[3] = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
					n = n.EuklNormalized();

					glm::vec3 from(0,0,1);
					glm::vec3 to(n[1],n[2],n[3]);
					q[ijk] = glm::quat(from,to);
				}
				else	// (0,0,0) is an invalid direction. Set to identity:
					q[ijk] = glm::quat(1,0,0,0);

				for(int d=0; d<stencil.nDir; d++)
					I[Index(ijk,d)] = stencil.W(d) * initialE[ijk];
			}
			else
			{// Kent intensity distribution:
			 // https://en.wikipedia.org/wiki/Kent_distribution
				// In case of dynamic streaming the stencil is rotated such that its north pole points towards n, whenever cxyz is used.
				// This means that the Kent distribution needs to point north for the dynamic streaming case.
			 	Tensor3 n = (isDynamicStreaming) ? Tensor3(0,0,1) : Tensor3(initialNx[ijk],initialNy[ijk],initialNz[ijk]);
				glm::vec3 from(0,0,1);
				glm::vec3 to(initialNx[ijk],initialNy[ijk],initialNz[ijk]);
				q[ijk] = glm::quat(from,to);

				for(size_t d1=0; d1<stencil.nPh; d1++)
				for(size_t d0=0; d0<stencil.nTh; d0++)
				{
					size_t d = stencil.Index(d0,d1);
					Tensor3 p = stencil.C(d0,d1);
					I[Index(ijk,d)] = stencil.W(d0,d1) * initialE[ijk] * exp(sigma * Tensor3::Dot(n, p));
				}
			}
		}
	}
}



void Radiation::NormalizeInitialIntensities()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
		if(isInitialGridPoint[ijk])
			initialE[ijk] *= initialE[ijk] / E_LF[ijk];
}



void Radiation::UpdateSphericalHarmonicsCoefficients()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	{
		size_t ijk = grid.Index(i,j,k);
		double dataS[lebedevStencil.nDir];
		double dataX[lebedevStencil.nDir];
		double dataY[lebedevStencil.nDir];
		double dataZ[lebedevStencil.nDir];
    	double dataCx[lebedevStencil.nDir];
    	double dataCy[lebedevStencil.nDir];
    	double dataCz[lebedevStencil.nDir];
		Coord xyz0 = grid.xyz(i,j,k);
		double alpha = metric.GetAlpha(ijk);

		for(size_t d=0; d<lebedevStencil.nDir; d++)
    	{
			// Initial data for geodesic equation:
			double s = 1;
			Coord xyz = xyz0;
            Tensor3 c = lebedevStencil.C(d);
            Tensor4 uIF(alpha, c[1] * alpha, c[2] * alpha, c[3] * alpha);
            Tensor3 vLF = Vec3ObservedByEulObs<IF,LF>(uIF, xyz, metric);

			// Solve geodesic equation backwards:
			if(!metric.InsideBH(xyz))
    	    	s *= RK45_GeodesicEquation<-1>(grid.dt, xyz, vLF, metric);
			else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
				vLF = Tensor3(0.0);

			Tensor3 vIF = TransformLFtoIF(vLF,metric.GetTetradInverse(xyz));

			// Final data points for fourier expansion:
			dataS[d] = 1.0/s;
			dataX[d] = xyz[1];
			dataY[d] = xyz[2];
			dataZ[d] = xyz[3];
			dataCx[d] = vIF[1];
			dataCy[d] = vIF[2];
			dataCz[d] = vIF[3];
    	}
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataS , &coefficientsS [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataX , &coefficientsX [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataY , &coefficientsY [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataZ , &coefficientsZ [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataCx, &coefficientsCx[HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataCy, &coefficientsCy[HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(lebedevStencil, dataCz, &coefficientsCz[HarmonicIndex(0,ijk)]);
	}
}



void Radiation::ComputeMomentsIF()
{
	PROFILE_FUNCTION();
	// constexpr double fourPiInv = 1.0 / (4.0 * M_PI);
	
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
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
		for(size_t d1=0; d1<stencil.nPh; d1++)
		for(size_t d0=0; d0<stencil.nTh; d0++)
		{
			// Tensor3 cxyz = stencil.C(d0,d1);
			Tensor3 cxyz = q[ijk] * stencil.C(d0,d1);
			size_t d = stencil.Index(d0,d1);
			size_t index = Index(ijk,d);
			E[ijk]   += I[index];
			Fx[ijk]  += I[index] * cxyz[1];
			Fy[ijk]  += I[index] * cxyz[2];
			Fz[ijk]  += I[index] * cxyz[3];
			Pxx[ijk] += I[index] * cxyz[1] * cxyz[1];
			Pxy[ijk] += I[index] * cxyz[1] * cxyz[2];
			Pxz[ijk] += I[index] * cxyz[1] * cxyz[3];
			Pyy[ijk] += I[index] * cxyz[2] * cxyz[2];
			Pyz[ijk] += I[index] * cxyz[2] * cxyz[3];
			Pzz[ijk] += I[index] * cxyz[3] * cxyz[3];
		}
        // E[ijk]   *= fourPiInv;
		// Tensor3 F = q[ijk] * Tensor3(Fx[ijk], Fy[ijk], Fz[ijk]);
        // Fx[ijk]  = F[1] * fourPiInv;
        // Fy[ijk]  = F[2] * fourPiInv;
        // Fz[ijk]  = F[3] * fourPiInv;
        // Fx[ijk]  *= fourPiInv;
        // Fy[ijk]  *= fourPiInv;
        // Fz[ijk]  *= fourPiInv;
        // Pxx[ijk] *= fourPiInv;
        // Pxy[ijk] *= fourPiInv;
        // Pxz[ijk] *= fourPiInv;
        // Pyy[ijk] *= fourPiInv;
        // Pyz[ijk] *= fourPiInv;
        // Pzz[ijk] *= fourPiInv;
	}
}
void Radiation::ComputeMomentsLF()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
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



Coord Radiation::GetTempCoordinate(size_t ijk, Tensor3 direction)
{
	Coord xyzTemp;
	xyzTemp[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsX[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
	xyzTemp[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsY[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
	xyzTemp[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsZ[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
	return xyzTemp;
}
Tensor3 Radiation::GetTemp3VelocityIF(size_t ijk, Tensor3 direction)
{
	Tensor3 vTempIF;
	vTempIF[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCx[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
	vTempIF[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCy[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
	vTempIF[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCz[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
	return vTempIF;
}
double Radiation::GetFrequencyShift(size_t ijk, Tensor3 direction)
{
	return SphericalHarmonicsXyz::GetValue(direction, &coefficientsS[HarmonicIndex(0,ijk)], lebedevStencil.nCoefficients);
}


// Linear:
//double Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF)
//{
//	// vTempIF is given in world space. Transform it to local space by applying inverse quaternion:
//	vTempIF = Invert(q[ijk]) * vTempIF;
//
//	// Index on local sphere grid of point (theta,phi).
//	double th = stencil.d0(vTempIF.Theta());
//	double ph = stencil.d1(vTempIF.Phi());
//
//	// Indices of nearest directions on sphere. Note that phi is cyclic and theta is not.
//	// These numbers are int instead of size_t due to possible negative values.
//	int th0 = floor(th);
//	int th1 = th0 + 1;
//	int ph0 = ((int)floor(ph) + stencil.nPh) % stencil.nPh;
//	int ph1 = (ph0 + 1) % stencil.nPh;
//
//	// Get intensities at nearest directions. Use north pole for th0==-1 and south pole for th1==nTh.
//	double I00 = (th0 == -1)          ? Inorth[ijk] : I[Index(ijk,th0,ph0)];
//	double I01 = (th0 == -1)          ? Inorth[ijk] : I[Index(ijk,th0,ph1)];
//	double I10 = (th1 == stencil.nTh) ? Isouth[ijk] : I[Index(ijk,th1,ph0)];
//	double I11 = (th1 == stencil.nTh) ? Isouth[ijk] : I[Index(ijk,th1,ph1)];
//
//	// Fractional part of sphere grid index.
//	// At poles it must be remapped to [0,1] because the index 0 and nTh-1
//	// do not correspond to north and south pole.
//	double tfrac = th - floor(th);
//	double pfrac = ph - floor(ph);
//	double d00 = stencil.d0(0);
//	if (th0 == -1)					// North
//		tfrac = 1.0 / abs(d00) * (tfrac - (1.0 + d00));
//	else if (th1 == stencil.nTh)	// South
//		tfrac /= abs(d00);
//
//	double value = BilinearInterpolation(tfrac, pfrac, I00, I01, I10, I11);
//	return std::max(value,0.0);
//}
//double Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF)
//{
//	// vTempIF is given in world space. Transform it to local space by applying inverse quaternion:
//	vTempIF = Invert(q[ijk]) * vTempIF;
//
//	// Index on local sphere grid of point (theta,phi).
//	double th = stencil.d0(vTempIF.Theta());
//	double ph = stencil.d1(vTempIF.Phi());
//
//	// Indices of nearest directions on sphere. Note that phi is cyclic and theta is not.
//	// These numbers are int instead of size_t due to possible negative values.
//	int th0 = floor(th);
//	int th1 = th0 + 1;
//	int ph0 = ((int)floor(ph) + stencil.nPh) % stencil.nPh;
//	int ph1 = (ph0 + 1) % stencil.nPh;
//
//    if (th0 == -1)
//    {
//        int d0 = th0 + 1;
//        double value = 0;
//        double sumInvDist = 0;
//        for (int d1=0; d1<stencil.nPh; d1++)
//        {
//            double invDist =  1.0 / Tensor3::UnitSphereNorm(vTempIF, stencil.C(d0,d1));
//            invDist = invDist * invDist;
//            value += I[Index(ijk,d0,d1)] * invDist;
//            sumInvDist += invDist;
//        }
//        value /= sumInvDist;
//        return std::max(value / stencil.W(d0,0), 0.0);
//    }
//    if (th0 == stencil.nTh)
//    {
//        int d0 = th0 - 1;
//        double value = 0;
//        double sumInvDist = 0;
//        for (int d1=0; d1<stencil.nPh; d1++)
//        {
//            double invDist =  1.0 / Tensor3::UnitSphereNorm(vTempIF, stencil.C(d0,d1));
//            invDist = invDist * invDist;
//            value += I[Index(ijk,d0,d1)] * invDist;
//            sumInvDist += invDist;
//        }
//        value /= sumInvDist;
//        return std::max(value / stencil.W(d0,0), 0.0);
//    }
//
//	// Interpolate Intensities from nearest stencil directions:
//	double I00 = I[Index(ijk,th0,ph0)] / stencil.W(th0,ph0);
//	double I01 = I[Index(ijk,th0,ph1)] / stencil.W(th0,ph1);
//	double I10 = I[Index(ijk,th1,ph0)] / stencil.W(th1,ph0);
//	double I11 = I[Index(ijk,th1,ph1)] / stencil.W(th1,ph1);
//	double value = BilinearInterpolation(th - floor(th), ph - floor(ph), I00, I01, I10, I11);
//	return std::max(value, 0.0);
//}
double Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF)
{
	// vTempIF is given in world space. Transform it to local space by applying inverse quaternion:
	vTempIF = Invert(q[ijk]) * vTempIF;

	// Index on local sphere grid of point (theta,phi).
	double th = stencil.d0(vTempIF.Theta());
	double ph = stencil.d1(vTempIF.Phi());
	int thP0 = floor(th);
	int phP0 = ((int)floor(ph) + stencil.nPh) % stencil.nPh;

	// Fractional part of sphere grid index.
	// At poles it must be remapped to [0,1] because the index 0 and nTh-1
	// do not correspond to north and south pole.
	double thFrac = th - floor(th);
	double phFrac = ph - floor(ph);

    double value;

    // North:
    if (thP0 == -1)
    {
        int thP1 = thP0 + 1;
        int phP1 = (phP0 + 1) % stencil.nPh;

        double I00 = Inorth[ijk];
        double I01 = Inorth[ijk];
        double I10 = I[Index(ijk,thP1,phP0)] / stencil.W(thP1,phP0);
        double I11 = I[Index(ijk,thP1,phP1)] / stencil.W(thP1,phP1);

    	double d00 = stencil.d0(0);
        thFrac = 1.0 / abs(d00) * (thFrac - (1.0 + d00));
        value = BilinearInterpolation(thFrac, phFrac, I00, I01, I10, I11);
    }
    // South:
    else if (thP0 == stencil.nTh - 1)
    {
        int thP1 = thP0 + 1;
        int phP1 = (phP0 + 1) % stencil.nPh;

        double I00 = I[Index(ijk,thP0,phP0)] / stencil.W(thP0,phP0);
        double I01 = I[Index(ijk,thP0,phP1)] / stencil.W(thP0,phP1);
        double I10 = Isouth[ijk];
        double I11 = Isouth[ijk];

    	double d00 = stencil.d0(0);
        thFrac /= abs(d00);
        value = BilinearInterpolation(thFrac, phFrac, I00, I01, I10, I11);
    }
    // North Band:
    else if (thP0 == 0)
    {
        int thP1 = thP0 + 1;
        int ipP2 = thP0 + 2;
        int phM1 = (phP0 - 1 + stencil.nPh) % stencil.nPh;
        int phP1 = (phP0 + 1) % stencil.nPh;
        int phP2 = (phP0 + 2) % stencil.nPh;

        double Im1 = SquaredInterpolation(thFrac - 1, I[Index(ijk,thP0,phM1)] / stencil.W(thP0,phM1), I[Index(ijk,thP1,phM1)] / stencil.W(thP1,phM1), I[Index(ijk,ipP2,phM1)] / stencil.W(ipP2,phM1));
        double Ip0 = SquaredInterpolation(thFrac - 1, I[Index(ijk,thP0,phP0)] / stencil.W(thP0,phP0), I[Index(ijk,thP1,phP0)] / stencil.W(thP1,phP0), I[Index(ijk,ipP2,phP0)] / stencil.W(ipP2,phP0));
        double Ip1 = SquaredInterpolation(thFrac - 1, I[Index(ijk,thP0,phP1)] / stencil.W(thP0,phP1), I[Index(ijk,thP1,phP1)] / stencil.W(thP1,phP1), I[Index(ijk,ipP2,phP1)] / stencil.W(ipP2,phP1));
        double Ip2 = SquaredInterpolation(thFrac - 1, I[Index(ijk,thP0,phP2)] / stencil.W(thP0,phP2), I[Index(ijk,thP1,phP2)] / stencil.W(thP1,phP2), I[Index(ijk,ipP2,phP2)] / stencil.W(ipP2,phP2));

        value = CubicInterpolation(phFrac, Im1, Ip0, Ip1, Ip2);
    }
    // South Band:
    else if (thP0 == stencil.nTh - 2)
    {
        int thM1 = thP0 - 1;
        int thP1 = thP0 + 1;
        int phM1 = (phP0 - 1 + stencil.nPh) % stencil.nPh;
        int phP1 = (phP0 + 1) % stencil.nPh;
        int phP2 = (phP0 + 2) % stencil.nPh;

        double Im1 = SquaredInterpolation(thFrac, I[Index(ijk,thM1,phM1)] / stencil.W(thM1,phM1), I[Index(ijk,thP0,phM1)] / stencil.W(thP0,phM1), I[Index(ijk,thP1,phM1)] / stencil.W(thP1,phM1));
        double Ip0 = SquaredInterpolation(thFrac, I[Index(ijk,thM1,phP0)] / stencil.W(thM1,phP0), I[Index(ijk,thP0,phP0)] / stencil.W(thP0,phP0), I[Index(ijk,thP1,phP0)] / stencil.W(thP1,phP0));
        double Ip1 = SquaredInterpolation(thFrac, I[Index(ijk,thM1,phP1)] / stencil.W(thM1,phP1), I[Index(ijk,thP0,phP1)] / stencil.W(thP0,phP1), I[Index(ijk,thP1,phP1)] / stencil.W(thP1,phP1));
        double Ip2 = SquaredInterpolation(thFrac, I[Index(ijk,thM1,phP2)] / stencil.W(thM1,phP2), I[Index(ijk,thP0,phP2)] / stencil.W(thP0,phP2), I[Index(ijk,thP1,phP2)] / stencil.W(thP1,phP2));

        value = CubicInterpolation(phFrac, Im1, Ip0, Ip1, Ip2);
    }
    // Bulk:
    else
    {
        int thM1 = thP0 - 1;
        int thP1 = thP0 + 1;
        int thP2 = thP0 + 2;
        int phM1 = (phP0 - 1 + stencil.nPh) % stencil.nPh;
        int phP1 = (phP0 + 1) % stencil.nPh;
        int phP2 = (phP0 + 2) % stencil.nPh;

        double Im1m1 = I[Index(ijk,thM1,phM1)] / stencil.W(thM1,phM1);
        double Im1p0 = I[Index(ijk,thM1,phP0)] / stencil.W(thM1,phP0);
        double Im1p1 = I[Index(ijk,thM1,phP1)] / stencil.W(thM1,phP1);
        double Im1p2 = I[Index(ijk,thM1,phP2)] / stencil.W(thM1,phP2);
        double Ip0m1 = I[Index(ijk,thP0,phM1)] / stencil.W(thP0,phM1);
        double Ip0p0 = I[Index(ijk,thP0,phP0)] / stencil.W(thP0,phP0);
        double Ip0p1 = I[Index(ijk,thP0,phP1)] / stencil.W(thP0,phP1);
        double Ip0p2 = I[Index(ijk,thP0,phP2)] / stencil.W(thP0,phP2);
        double Ip1m1 = I[Index(ijk,thP1,phM1)] / stencil.W(thP1,phM1);
        double Ip1p0 = I[Index(ijk,thP1,phP0)] / stencil.W(thP1,phP0);
        double Ip1p1 = I[Index(ijk,thP1,phP1)] / stencil.W(thP1,phP1);
        double Ip1p2 = I[Index(ijk,thP1,phP2)] / stencil.W(thP1,phP2);
        double Ip2m1 = I[Index(ijk,thP2,phM1)] / stencil.W(thP2,phM1);
        double Ip2p0 = I[Index(ijk,thP2,phP0)] / stencil.W(thP2,phP0);
        double Ip2p1 = I[Index(ijk,thP2,phP1)] / stencil.W(thP2,phP1);
        double Ip2p2 = I[Index(ijk,thP2,phP2)] / stencil.W(thP2,phP2);
        value = BicubicInterpolation
                (thFrac, phFrac,
                 Im1m1, Im1p0, Im1p1, Im1p2,
                 Ip0m1, Ip0p0, Ip0p1, Ip0p2,
                 Ip1m1, Ip1p0, Ip1p1, Ip1p2,
                 Ip2m1, Ip2p0, Ip2p1, Ip2p2);
    }
    return std::max(value, 0.0);
}



Tensor3 Radiation::AverageF(size_t i, size_t j, size_t k)
{
	// i,j,k >= 2, thus i+a etc will never be negative.
	Tensor3 averageF(0.0);
	for(int c=-1; c<=1; c++)
	for(int b=-1; b<=1; b++)
	for(int a=-1; a<=1; a++)
	{
		size_t index = grid.Index(i+a, j+b, k+c);
		averageF[1] += Fx[index];
		averageF[2] += Fy[index];
		averageF[3] += Fz[index];
	}
	averageF[1] /= 27.0;
	averageF[2] /= 27.0;
	averageF[3] /= 27.0;
	return averageF;
}



void Radiation::SetPoleIntensities()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
    {
        int ijk = grid.Index(i,j,k);

        Inorth[ijk] = Isouth[ijk] = 0;
        int d0north = 0;
        int d0south = stencil.nTh - 1;
        for(int d1=0; d1<stencil.nPh; d1++)
        {
            Inorth[ijk] += I[Index(ijk,d0north,d1)] / stencil.W(d0north,d1);
            Isouth[ijk] += I[Index(ijk,d0south,d1)] / stencil.W(d0south,d1);
        }
        Inorth[ijk] /= stencil.nPh;
        Isouth[ijk] /= stencil.nPh;
    }
}



void Radiation::UpdateQuaternions()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	{
		size_t ijk = grid.Index(i,j,k);
		Tensor3 averageF = AverageF(i,j,k);
		double norm = averageF.EuklNorm();

		// At least 1% of the lights points in the direction of the first momentum
		if (norm / E[ijk] > 0.01)
		{
			glm::vec3 from(0,0,1);
			glm::vec3 to(averageF[1]/norm,averageF[2]/norm,averageF[3]/norm);
			qNew[ijk] = glm::quat(from,to);
		}
		else
			qNew[ijk] = q[ijk];
	}
}



void Radiation::StreamFlatStatic()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(5)
	#ifdef ijkd
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	#endif
	#ifdef dijk
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	#endif
		StreamFlatKernal<Static>(i,j,k,d0,d1);

	std::swap(I,Inew);
}
void Radiation::StreamFlatDynamic()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(5)
	#ifdef ijkd
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	#endif
	#ifdef dijk
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	#endif
		StreamFlatKernal<Dynamic>(i,j,k,d0,d1);

	std::swap(I,Inew);
	std::swap(q,qNew);
}
template<class StaticOrDynamic>
void Radiation::StreamFlatKernal(size_t i, size_t j,size_t k, size_t d0, size_t d1)
{
	size_t ijk = grid.Index(i,j,k);	    // Index of lattice point ijk
	size_t d = stencil.Index(d0,d1);	// Index of direction d
	size_t index = Index(ijk,d);		// Index of population d at lattice point ijk

	// Get temp velocity:
	Tensor3 vTempIF;
	if constexpr(std::is_same<StaticOrDynamic,Dynamic>::value)
        vTempIF = qNew[ijk] * stencil.C(d0,d1);
	if constexpr(std::is_same<StaticOrDynamic,Static>::value)
		vTempIF = stencil.C(d0,d1);

	// Get temp lattice point:
	Coord xyzTemp = grid.xyz(i,j,k);
	xyzTemp[1] -= vTempIF[1] * grid.dt;
	xyzTemp[2] -= vTempIF[2] * grid.dt;
	xyzTemp[3] -= vTempIF[3] * grid.dt;

	// Get 8 nearest Grid Points:
	double iTemp = grid.i(xyzTemp[1]);
	double jTemp = grid.j(xyzTemp[2]);
	double kTemp = grid.k(xyzTemp[3]);
	size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
	size_t k0 = std::floor(kTemp);	size_t k1 = k0 + 1;

	// Interpolate intensity from neighbouring 8 lattice points to temporary point:
	if constexpr(std::is_same<StaticOrDynamic,Dynamic>::value)
	{
		double intensityAt_i0j0k0 = IntensityAt(grid.Index(i0,j0,k0),vTempIF);
		double intensityAt_i0j0k1 = IntensityAt(grid.Index(i0,j0,k1),vTempIF);
		double intensityAt_i0j1k0 = IntensityAt(grid.Index(i0,j1,k0),vTempIF);
		double intensityAt_i0j1k1 = IntensityAt(grid.Index(i0,j1,k1),vTempIF);
		double intensityAt_i1j0k0 = IntensityAt(grid.Index(i1,j0,k0),vTempIF);
		double intensityAt_i1j0k1 = IntensityAt(grid.Index(i1,j0,k1),vTempIF);
		double intensityAt_i1j1k0 = IntensityAt(grid.Index(i1,j1,k0),vTempIF);
		double intensityAt_i1j1k1 = IntensityAt(grid.Index(i1,j1,k1),vTempIF);
		Inew[index] = stencil.W(d0,d1) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
					  intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
					  intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
	}
	if constexpr(std::is_same<StaticOrDynamic,Static>::value)
		Inew[index] = TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		              I[Index(i0,j0,k0,d)], I[Index(i0,j0,k1,d)], I[Index(i0,j1,k0,d)], I[Index(i0,j1,k1,d)],
		              I[Index(i1,j0,k0,d)], I[Index(i1,j0,k1,d)], I[Index(i1,j1,k0,d)], I[Index(i1,j1,k1,d)]);
}



void Radiation::StreamCurvedStatic()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(5)
	#ifdef ijkd
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	#endif
	#ifdef dijk
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	#endif
		StreamCurvedKernal<Static>(i,j,k,d0,d1);

	std::swap(I,Inew);
}
void Radiation::StreamCurvedDynamic()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(5)
	#ifdef ijkd
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	#endif
	#ifdef dijk
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	for(size_t d1=0; d1<stencil.nPh; d1++)
	for(size_t d0=0; d0<stencil.nTh; d0++)
	#endif
		StreamCurvedKernal<Dynamic>(i,j,k,d0,d1);

	std::swap(I,Inew);
	std::swap(q,qNew);
}
template<class StaticOrDynamic>
void Radiation::StreamCurvedKernal(size_t i, size_t j, size_t k, size_t d0, size_t d1)
{
	size_t ijk = grid.Index(i,j,k);		// Index of lattice point ijk
	size_t d = stencil.Index(d0,d1);	// Index of direction d
	size_t index = Index(ijk,d);		// Index of population d at lattice point ijk

	// Skip LPs which are inside BH:
	if(metric.InsideBH(grid.xyz(i,j,k)))
	{
        Inew[index] = 0;
		return;
	}

	// Get velocity direction in IF:
	Tensor3 direction;
	if constexpr(std::is_same<StaticOrDynamic,Dynamic>::value)
        direction = qNew[ijk] * stencil.C(d0,d1);
	if constexpr(std::is_same<StaticOrDynamic,Static>::value)
        direction = stencil.C(d0,d1);
	
	// Get quantities at emission point:
	double s = GetFrequencyShift(ijk, direction);
	Coord xyzTemp = GetTempCoordinate(ijk, direction);
	Tensor3 vTempIF = GetTemp3VelocityIF(ijk, direction);

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
	size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
	size_t k0 = std::floor(kTemp);	size_t k1 = k0 + 1;

	// Intensity interpolation:
	double alpha = metric.GetAlpha(ijk);
	double intensityAt_i0j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k0))) * IntensityAt(grid.Index(i0,j0,k0),vTempIF);
	double intensityAt_i0j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k1))) * IntensityAt(grid.Index(i0,j0,k1),vTempIF);
	double intensityAt_i0j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k0))) * IntensityAt(grid.Index(i0,j1,k0),vTempIF);
	double intensityAt_i0j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k1))) * IntensityAt(grid.Index(i0,j1,k1),vTempIF);
	double intensityAt_i1j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k0))) * IntensityAt(grid.Index(i1,j0,k0),vTempIF);
	double intensityAt_i1j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k1))) * IntensityAt(grid.Index(i1,j0,k1),vTempIF);
	double intensityAt_i1j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k0))) * IntensityAt(grid.Index(i1,j1,k0),vTempIF);
	double intensityAt_i1j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k1))) * IntensityAt(grid.Index(i1,j1,k1),vTempIF);

	// Interpolate intensity from neighbouring 4 lattice points to temporary point:
	Inew[index] = stencil.W(d0,d1) * IntegerPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
				  intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
				  intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
}



void Radiation::Collide()
{
	PROFILE_FUNCTION();
	// TODO: Steife DGL?
	
	PARALLEL_FOR(3)
	for(size_t k = 2; k < metric.grid.nz - 2; k++)
	for(size_t j = 2; j < metric.grid.ny - 2; j++)
	for(size_t i = 2; i < metric.grid.nx - 2; i++)
	{
		if(metric.InsideBH(grid.xyz(i,j,k)))
			continue;

		size_t ijk = grid.Index(i,j,k);

		// Simulate stationary fluid, u^k=(0,0):
		double alpha = metric.GetAlpha(ijk);
		Tensor3 u(0.0);	// fluid 3 velocity as seen by Eulerian observer
		double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ijk)));	// Lorentz factor
		double uDotF = Fx_LF[ijk] * u[1] + Fy_LF[ijk] * u[2] + Fz_LF[ijk] * u[3];	// F^i u_i
		double uuDotP =
		Pxx_LF[ijk] * u[1] * u[1] + Pyy_LF[ijk] * u[2] * u[2] + Pzz_LF[ijk] * u[3] * u[3]
		+ 2.0 * (Pxy_LF[ijk] * u[1] * u[2] + Pxz_LF[ijk] * u[1] * u[3] + Pyz_LF[ijk] * u[2] * u[3]);	// P^ij u_i u_j
		double fluidE = W * W * (E_LF[ijk] - 2.0 * uDotF + uuDotP);
		
		for(size_t d1 = 0; d1 < stencil.nPh; d1++)
		for(size_t d0 = 0; d0 < stencil.nTh; d0++)
		{
			size_t d = stencil.Index(d0,d1);
			size_t index = Index(ijk,d);
			double A = W * (1.0 - Tensor3::Dot(q[ijk] * stencil.C(d0,d1), u));

			double Gamma = stencil.W(d0,d1) * (eta[ijk] + kappa0[ijk]*fluidE) / (A*A*A) - A*I[index] * (kappaA[ijk] + kappa0[ijk]);
			I[index] += alpha * metric.grid.dt * Gamma;
		}
	}
}



void Radiation::TakePicture()
{
	// Note:
	// ij is camera index
	// i,j,k are grid indexes and not related to ij.
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ij=0; ij<camera.pixelCount; ij++)
	{
		Coord pixel = camera.xyz(ij);
		if(grid.OutsideDomain(pixel))
		{
			camera.image[ij] = 0;
			continue;
		}
		
		// We want to measure only the intensities that move orthogonal through the camera plane,
		// meaning the light velocity is the opposite of the camera normal vector.
		Tensor3 lookDir = camera.lookDirection;
		Tensor4 uLF(1,-lookDir[1],-lookDir[2],-lookDir[3]);
		uLF = NullNormalize(uLF,metric.GetMetric_ll(pixel));
		Tensor3 vIF = Vec3ObservedByEulObs<LF,IF>(uLF, pixel, metric);

		// Get 8 nearest Grid Points:
		double iTemp = grid.i(pixel[1]);
		double jTemp = grid.j(pixel[2]);
		double kTemp = grid.k(pixel[3]);
		size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
		size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
		size_t k0 = std::floor(kTemp);	size_t k1 = k0 + 1;
		
		// Intensity interpolation:
		//double alpha = metric.GetAlpha(pixel); // removed to get intensity as seen by observer infinitly far away.
		double intensityAt_i0j0k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0,j0,k0))) * IntensityAt(grid.Index(i0,j0,k0),vIF);
		double intensityAt_i0j0k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0,j0,k1))) * IntensityAt(grid.Index(i0,j0,k1),vIF);
		double intensityAt_i0j1k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0,j1,k0))) * IntensityAt(grid.Index(i0,j1,k0),vIF);
		double intensityAt_i0j1k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i0,j1,k1))) * IntensityAt(grid.Index(i0,j1,k1),vIF);
		double intensityAt_i1j0k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1,j0,k0))) * IntensityAt(grid.Index(i1,j0,k0),vIF);
		double intensityAt_i1j0k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1,j0,k1))) * IntensityAt(grid.Index(i1,j0,k1),vIF);
		double intensityAt_i1j1k0 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1,j1,k0))) * IntensityAt(grid.Index(i1,j1,k0),vIF);
		double intensityAt_i1j1k1 = IntegerPow<4>(1.0 / metric.GetAlpha(grid.Index(i1,j1,k1))) * IntensityAt(grid.Index(i1,j1,k1),vIF);	
		
		camera.image[ij]
		= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
		 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);

		// Energy density interpolation:
		//camera.image[ij]
		//= TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
		// E_LF[grid.Index(i0,j0,k0)], E_LF[grid.Index(i0,j0,k1)], E_LF[grid.Index(i0,j1,k0)], E_LF[grid.Index(i0,j1,k1)],
		// E_LF[grid.Index(i1,j0,k0)], E_LF[grid.Index(i1,j0,k1)], E_LF[grid.Index(i1,j1,k0)], E_LF[grid.Index(i1,j1,k1)]);
	}
}



void Radiation::WriteIntensitiesToCsv(double time, const int frameNumber, std::string directory, std::string name)
{
	PROFILE_FUNCTION();
    CreateDirectory(directory);

    name = name + name + FrameNumber(frameNumber) + ".csv";
    std::ofstream fileOut(directory + "/" + name);

    double Imin =  1e20;
    double Imax = -1e20;
	for(size_t k = 2; k < metric.grid.nz - 2; k++)
	for(size_t j = 2; j < metric.grid.ny - 2; j++)
	for(size_t i = 2; i < metric.grid.nx - 2; i++)
    for(size_t d=0; d<stencil.nDir; d++)
    {
        int index = Index(i,j,k,d);
        Imin = std::min(Imin, I[index]);
        Imax = std::max(Imax, I[index]);
    }


	fileOut << "#x, y, z, color\n";
	for(size_t k = 2; k < metric.grid.nz - 2; k++)
	for(size_t j = 2; j < metric.grid.ny - 2; j++)
	for(size_t i = 2; i < metric.grid.nx - 2; i++)
	{
		size_t ijk = grid.Index(i,j,k);
		Coord xyz = grid.xyz(i,j,k);
		fileOut << xyz[1] << ", " << xyz[2] << ", " << xyz[3] << ", " << 0 << "\n";
		// Bulk:
		for(size_t d=0; d<stencil.nDir; d++)
		{
			size_t index = Index(ijk,d);
			if(I[index] > 1e-8)
			{
				Tensor3 dir = q[ijk] * stencil.C(d);
                double value = (I[index] - Imin) / (Imax - Imin) / stencil.W(d);

				Coord pos = xyz;
				pos[1] += dir[1] * grid.dt * 0.5;// * (0.1 + 0.9 * value);
				pos[2] += dir[2] * grid.dt * 0.5;// * (0.1 + 0.9 * value);
				pos[3] += dir[3] * grid.dt * 0.5;// * (0.1 + 0.9 * value);

				fileOut << pos[1] << ", " << pos[2] << ", " << pos[3] << ", " << value << "\n";
			}
		}
		// Orthogonal to Camera:
		{
			Tensor3 lookDir = camera.lookDirection;
			Tensor4 uLF(1,-lookDir[1],-lookDir[2],-lookDir[3]);
			uLF = NullNormalize(uLF,metric.GetMetric_ll(xyz));
			Tensor3 dir = Vec3ObservedByEulObs<LF,IF>(uLF, xyz, metric);
			
			double I = IntensityAt(ijk,dir);
            double value = (I - Imin) / (Imax - Imin);

			Coord pos = xyz;
			pos[1] += dir[1] * grid.dt * 0.7;// * (0.1 + 0.9 * value);
			pos[2] += dir[2] * grid.dt * 0.7;// * (0.1 + 0.9 * value);
			pos[3] += dir[3] * grid.dt * 0.7;// * (0.1 + 0.9 * value);
			
			fileOut << pos[1] << ", " << pos[2] << ", " << pos[3] << ", " << value << "\n";
		}
	}
	
    fileOut.close();
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
		std::cout << " sigma        = " << sigma << "\n";
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

			if((config.writeData || config.useCamera) && (n % config.writeFrequency) == 0)
				ComputeMomentsLF();
			if(n % config.writeFrequency == 0)
			{
				if(config.writeData)
				{
					grid.WriteFrametoCsv(n*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, n, logger.directoryPath + "/Moments", config.name);
					WriteIntensitiesToCsv(n*grid.dt, n, logger.directoryPath + "/Intensities", "I");
				}
				if(config.useCamera && (n % config.writeFrequency) == 0)
				{
					TakePicture();
					camera.WriteImagetoCsv(n*grid.dt, n, logger.directoryPath + "/CameraImages");
				}
			}

			if (config.updateSphericalHarmonics)
			{ UpdateSphericalHarmonicsCoefficients(); }

			// Streaming:
            SetPoleIntensities();
			switch(streamingType)
			{
				case(StreamingType::FlatStatic):		StreamFlatStatic();							break;
				case(StreamingType::FlatDynamic):		UpdateQuaternions(); StreamFlatDynamic();	break;
				case(StreamingType::CurvedStatic):		StreamCurvedStatic();	 					break;
				case(StreamingType::CurvedDynamic):		UpdateQuaternions(); StreamCurvedDynamic();	break;
			}

			if(config.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if(config.writeData || config.useCamera)
			ComputeMomentsLF();
		if(config.writeData)
			grid.WriteFrametoCsv(timeSteps*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, timeSteps, logger.directoryPath + "/Moments", config.name);
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
*/