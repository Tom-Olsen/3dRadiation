#include "Radiation.h"

#define HALO 1

Radiation::Radiation(Metric& metric, Stencil& stencil, LebedevStencil& streamingStencil, Camera& camera, StreamingType streamingType):
grid(metric.grid), metric(metric), stencil(stencil), streamingStencil(streamingStencil), camera(camera), streamingType(streamingType)
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
	coefficientsS.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsX.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsY.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsZ.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsCx.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsCy.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsCz.resize(grid.nxyz * streamingStencil.nCoefficients);
    
	// Initialize all Quaternions to identity:
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
        q[ijk] = qNew[ijk] = glm::quat(1,0,0,0);
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
		return d + ijk * stencil.nDir;
	#endif
}
size_t Radiation::Index(size_t i, size_t j, size_t k, size_t d)
{
	size_t ijk = grid.Index(i,j,k);
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return d + ijk * stencil.nDir;
	#endif
}
size_t Radiation::HarmonicIndex(size_t f, size_t ijk)
{
	return f + ijk * streamingStencil.nCoefficients;
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
	bool isAdaptiveStreaming = (streamingType == StreamingType::FlatAdaptive || streamingType == StreamingType::CurvedAdaptive);
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
				if (isAdaptiveStreaming)
				{// Random direction for uniformity:
					Tensor3 n;
                    double norm = 2;
                    while(norm > 1)
                    {
					    n[1] = RandomRange(-1,1);
					    n[2] = RandomRange(-1,1);
					    n[3] = RandomRange(-1,1);
                        norm = n.EuklNorm();
                    }
                
					glm::vec3 to(n[1]/norm, n[2]/norm, n[3]/norm);
					q[ijk] = glm::quat(from,to);
				}
				else	// (0,0,0) is an invalid direction. Set to identity:
					q[ijk] = glm::quat(1,0,0,0);

				for(int d=0; d<stencil.nDir; d++)
					I[Index(ijk,d)] = initialE[ijk];
			}
			else
			{// Kent intensity distribution:
			 // https://en.wikipedia.org/wiki/Kent_distribution
				// In case of adaptive streaming the stencil is rotated such that its north pole points towards n, whenever cxyz is used.
				// This means that the Kent distribution needs to point north for the adaptive streaming case.
			 	Tensor3 n = (isAdaptiveStreaming) ? Tensor3(0,0,1) : Tensor3(initialNx[ijk],initialNy[ijk],initialNz[ijk]);
				glm::vec3 to(initialNx[ijk],initialNy[ijk],initialNz[ijk]);
				q[ijk] = glm::quat(from,to);

				for(size_t d=0; d<stencil.nDir; d++)
				{
					Tensor3 p = stencil.Ct3(d);
					I[Index(ijk,d)] = initialE[ijk] * exp(sigma * Tensor3::Dot(n, p));
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
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	{
		size_t ijk = grid.Index(i,j,k);
		double dataS[streamingStencil.nDir];
		double dataX[streamingStencil.nDir];
		double dataY[streamingStencil.nDir];
		double dataZ[streamingStencil.nDir];
    	double dataCx[streamingStencil.nDir];
    	double dataCy[streamingStencil.nDir];
    	double dataCz[streamingStencil.nDir];
		Coord xyz0 = grid.xyz(i,j,k);
		double alpha = metric.GetAlpha(ijk);

		for(size_t d=0; d<streamingStencil.nDir; d++)
    	{
			// Initial data for geodesic equation:
			double s = 1;
			Coord xyz = xyz0;
            Tensor3 c = streamingStencil.Ct3(d);
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
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataS , &coefficientsS [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataX , &coefficientsX [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataY , &coefficientsY [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataZ , &coefficientsZ [HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[HarmonicIndex(0,ijk)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCz, &coefficientsCz[HarmonicIndex(0,ijk)]);
	}
}



void Radiation::ComputeMomentsIF()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
	{
		E[ijk]   = 0.0;
		Fx[ijk]  = 0.0;
		Fy[ijk]  = 0.0;
		Fz[ijk]  = 0.0;
		Pxx[ijk] = 0.0;
		Pxy[ijk] = 0.0;
		Pxz[ijk] = 0.0;
		Pyy[ijk] = 0.0;
		Pyz[ijk] = 0.0;
		Pzz[ijk] = 0.0;
		for(size_t d=0; d<stencil.nDir; d++)
		{
			Tensor3 dir = q[ijk] * stencil.Ct3(d);
			size_t index = Index(ijk,d);
            double c = stencil.W(d) * I[index];
			E[ijk]   += c;
			Fx[ijk]  += c * dir[1];
			Fy[ijk]  += c * dir[2];
			Fz[ijk]  += c * dir[3];
			Pxx[ijk] += c * dir[1] * dir[1];
			Pxy[ijk] += c * dir[1] * dir[2];
			Pxz[ijk] += c * dir[1] * dir[3];
			Pyy[ijk] += c * dir[2] * dir[2];
			Pyz[ijk] += c * dir[2] * dir[3];
			Pzz[ijk] += c * dir[3] * dir[3];
		}
	}
}
void Radiation::ComputeMomentsLF()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=0; k<grid.nz; k++)
	for(size_t j=0; j<grid.ny; j++)
	for(size_t i=0; i<grid.nx; i++)
	{
        int ijk = grid.Index(i,j,k);
        if(metric.InsideBH(grid.xyz(i,j,k)))
        {
		    E_LF[ijk]   = 0.0;
		    Fx_LF[ijk]  = 0.0;
		    Fy_LF[ijk]  = 0.0;
		    Fz_LF[ijk]  = 0.0;
		    Pxx_LF[ijk] = 0.0;
		    Pxy_LF[ijk] = 0.0;
		    Pxz_LF[ijk] = 0.0;
		    Pyy_LF[ijk] = 0.0;
		    Pyz_LF[ijk] = 0.0;
		    Pzz_LF[ijk] = 0.0;
            continue;
        }
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



Coord Radiation::GetTempCoordinate(size_t ijk, const Tensor3& direction)
{
	Coord xyzTemp;
	xyzTemp[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsX[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
	xyzTemp[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsY[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
	xyzTemp[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsZ[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
	return xyzTemp;
}
Tensor3 Radiation::GetTemp3VelocityIF(size_t ijk, const Tensor3& direction)
{
	Tensor3 vTempIF;
	vTempIF[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCx[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
	vTempIF[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCy[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
	vTempIF[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCz[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
	return vTempIF;
}
double Radiation::GetFrequencyShift(size_t ijk, const Tensor3& direction)
{
	return SphericalHarmonicsXyz::GetValue(direction, &coefficientsS[HarmonicIndex(0,ijk)], streamingStencil.nCoefficients);
}



double Radiation::IntensityAt(size_t ijk, Tensor3 vTempIF)
{
    // Barycentric interpolation:
    // vTempIF = Invert(q[ijk]) * vTempIF;
    // Vector3Int triangle;
    // Vector3 weights = stencil.BarycentricWeights(vTempIF, triangle);
    // double I0 = I[Index(ijk,triangle[0])];
    // double I1 = I[Index(ijk,triangle[1])];
    // double I2 = I[Index(ijk,triangle[2])];
    // return std::max(0.0, I0 * weights[0] + I1 * weights[1] + I2 * weights[2]);

    // Inverse Distance Interpolation:
    // vTempIF = Invert(q[ijk]) * vTempIF;
    // size_t d = stencil.NearestNeighbour(vTempIF);
    // std::span naturalNeighbours = stencil.VoronoiNeighbourOf(d);
    // double value = 0;
    // double invDistSum = 0;
    // for(int k : naturalNeighbours)
    // {
        // double dist = Tensor3::UnitSphereNorm(vTempIF,stencil.Ct3(k));
        // if(dist < 1e-16)    // vTempIF == stencil.Ct3(d) == stencil.Ct3(k)
            // return I[Index(ijk,k)];
        // double invDist = 1.0 / dist;
        // value += I[Index(ijk,k)] * invDist;
        // invDistSum += invDist;
    // }
    // return std::max(0.0, value/invDistSum);

    // Voronoi Interpolation:
    // vTempIF = Invert(q[ijk]) * vTempIF;
    // std::vector<int> neighbours;
    // std::vector<double> weights = stencil.VoronoiWeights(vTempIF, neighbours);
    // double value = 0.0;
    // for(int i=0; i<weights.size(); i++)
        // value += weights[i] * I[Index(ijk,neighbours[i])];
    // return std::max(0.0, value);
    
    // Voronoi Interpolation, optimized:
    vTempIF = Invert(q[ijk]) * vTempIF;
    double i = stencil.interpolationGrid.i(vTempIF.Theta());
    double j = std::fmod((stencil.interpolationGrid.j(vTempIF.Phi()) + stencil.interpolationGrid.nPh), stencil.interpolationGrid.nPh);
    size_t i0 = std::floor(i); size_t i1 = i0 + 1;
    size_t j0 = std::floor(j); size_t j1 = (j0 + 1) % stencil.interpolationGrid.nPh;
    std::span<const size_t> neighbours00 = stencil.voronoiNeighboursOnGrid.Row(stencil.interpolationGrid.Index(i0,j0));
    std::span<const size_t> neighbours01 = stencil.voronoiNeighboursOnGrid.Row(stencil.interpolationGrid.Index(i0,j1));
    std::span<const size_t> neighbours10 = stencil.voronoiNeighboursOnGrid.Row(stencil.interpolationGrid.Index(i1,j0));
    std::span<const size_t> neighbours11 = stencil.voronoiNeighboursOnGrid.Row(stencil.interpolationGrid.Index(i1,j1));
    std::span<const double> weights00 = stencil.voronoiWeightsOnGrid.Row(stencil.interpolationGrid.Index(i0,j0));
    std::span<const double> weights01 = stencil.voronoiWeightsOnGrid.Row(stencil.interpolationGrid.Index(i0,j1));
    std::span<const double> weights10 = stencil.voronoiWeightsOnGrid.Row(stencil.interpolationGrid.Index(i1,j0));
    std::span<const double> weights11 = stencil.voronoiWeightsOnGrid.Row(stencil.interpolationGrid.Index(i1,j1));
    double value00 = 0;
    double value01 = 0;
    double value10 = 0;
    double value11 = 0;
    for(size_t i=0; i<weights00.size(); i++)
        value00 += weights00[i] * I[Index(ijk,neighbours00[i])];
    for(size_t i=0; i<weights01.size(); i++)
        value01 += weights01[i] * I[Index(ijk,neighbours01[i])];
    for(size_t i=0; i<weights10.size(); i++)
        value10 += weights10[i] * I[Index(ijk,neighbours10[i])];
    for(size_t i=0; i<weights11.size(); i++)
        value11 += weights11[i] * I[Index(ijk,neighbours11[i])];
    return std::max(0.0, BilinearInterpolation(i-i0,j-j0,value00,value01,value10,value11));
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



void Radiation::UpdateQuaternions()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	{
		size_t ijk = grid.Index(i,j,k);
		Tensor3 averageF = AverageF(i,j,k);
		double norm = averageF.EuklNorm();

		// At least 1% of the light points in the direction of the first momentum
		if (E[ijk] > 1e-16 && norm / E[ijk] > 0.01)
		{
			glm::vec3 to(averageF[1]/norm, averageF[2]/norm, averageF[3]/norm);
			qNew[ijk] = glm::quat(from,to);
		}
		else
			qNew[ijk] = q[ijk];
	}
}



void Radiation::StreamFlatFixed()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(4)
	#ifdef ijkd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dijk
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
	    size_t ijk = grid.Index(i,j,k);	    // Index of lattice point ijk
	    size_t index = Index(ijk,d);		// Index of population d at lattice point ijk

	    // Get temp velocity:
	    Tensor3 direction = stencil.Ct3(d);

	    // Get temp lattice point:
	    Coord xyzTemp = grid.xyz(i,j,k);
	    xyzTemp[1] -= direction[1] * grid.dt;
	    xyzTemp[2] -= direction[2] * grid.dt;
	    xyzTemp[3] -= direction[3] * grid.dt;

	    // Get 8 nearest Grid Points:
	    double iTemp = grid.i(xyzTemp[1]);
	    double jTemp = grid.j(xyzTemp[2]);
	    double kTemp = grid.k(xyzTemp[3]);
	    size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	    size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
	    size_t k0 = std::floor(kTemp);	size_t k1 = k0 + 1;
        
	    Inew[index] = TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
	                I[Index(i0,j0,k0,d)], I[Index(i0,j0,k1,d)], I[Index(i0,j1,k0,d)], I[Index(i0,j1,k1,d)],
	                I[Index(i1,j0,k0,d)], I[Index(i1,j0,k1,d)], I[Index(i1,j1,k0,d)], I[Index(i1,j1,k1,d)]);
    }
	std::swap(I,Inew);
}

void Radiation::StreamFlatAdaptive()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(4)
	#ifdef ijkd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dijk
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
	    size_t ijk = grid.Index(i,j,k);	    // Index of lattice point ijk
	    size_t index = Index(ijk,d);		// Index of population d at lattice point ijk

	    // Get temp velocity:
	    Tensor3 direction = qNew[ijk] * stencil.Ct3(d);

	    // Get temp lattice point:
	    Coord xyzTemp = grid.xyz(i,j,k);
	    xyzTemp[1] -= direction[1] * grid.dt;
	    xyzTemp[2] -= direction[2] * grid.dt;
	    xyzTemp[3] -= direction[3] * grid.dt;

	    // Get 8 nearest Grid Points:
	    double iTemp = grid.i(xyzTemp[1]);
	    double jTemp = grid.j(xyzTemp[2]);
	    double kTemp = grid.k(xyzTemp[3]);
	    size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	    size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
	    size_t k0 = std::floor(kTemp);	size_t k1 = k0 + 1;

		double intensityAt_i0j0k0 = IntensityAt(grid.Index(i0,j0,k0),direction);
		double intensityAt_i0j0k1 = IntensityAt(grid.Index(i0,j0,k1),direction);
		double intensityAt_i0j1k0 = IntensityAt(grid.Index(i0,j1,k0),direction);
		double intensityAt_i0j1k1 = IntensityAt(grid.Index(i0,j1,k1),direction);
		double intensityAt_i1j0k0 = IntensityAt(grid.Index(i1,j0,k0),direction);
		double intensityAt_i1j0k1 = IntensityAt(grid.Index(i1,j0,k1),direction);
		double intensityAt_i1j1k0 = IntensityAt(grid.Index(i1,j1,k0),direction);
		double intensityAt_i1j1k1 = IntensityAt(grid.Index(i1,j1,k1),direction);
		Inew[index] = TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
					  intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
					  intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
    }
	std::swap(I,Inew);
	std::swap(q,qNew);
}



void Radiation::StreamCurvedFixed()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(4)
	#ifdef ijkd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dijk
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
        size_t ijk = grid.Index(i,j,k);		// Index of lattice point ijk
	    size_t index = Index(ijk,d);		// Index of population d at lattice point ijk

	    // Skip LPs which are inside BH:
	    if(metric.InsideBH(grid.xyz(i,j,k)))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get velocity direction in IF:
	    Tensor3 direction = stencil.Ct3(d);
    
	    // Get quantities at emission point:
	    double s = GetFrequencyShift(ijk, direction);
	    Coord xyzTemp = GetTempCoordinate(ijk, direction);
	    Tensor3 vTempIF = GetTemp3VelocityIF(ijk, direction);

	    // Skip temporary Grid Points inside BH:
	    if(metric.InsideBH(xyzTemp))
	    {
            Inew[index] = 0;
	    	continue;
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

	    // Interpolate intensity from neighbouring 8 lattice points to temporary point:
	    Inew[index] = IntegerPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
	    			  intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
	    			  intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
    }
	std::swap(I,Inew);
}

void Radiation::StreamCurvedAdaptive()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(4)
	#ifdef ijkd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dijk
	for(size_t k=HALO; k<grid.nz-HALO; k++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
        size_t ijk = grid.Index(i,j,k);		// Index of lattice point ijk
	    size_t index = Index(ijk,d);		// Index of population d at lattice point ijk

	    // Skip LPs which are inside BH:
	    if(metric.InsideBH(grid.xyz(i,j,k)))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get velocity direction in IF:
	    Tensor3 direction = qNew[ijk] * stencil.Ct3(d);
    
	    // Get quantities at emission point:
	    double s = GetFrequencyShift(ijk, direction);
	    Coord xyzTemp = GetTempCoordinate(ijk, direction);
	    Tensor3 vTempIF = GetTemp3VelocityIF(ijk, direction);

	    // Skip temporary Grid Points inside BH:
	    if(metric.InsideBH(xyzTemp))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get 8 nearest Grid Points:
	    double iTemp = grid.i(xyzTemp[1]);
	    double jTemp = grid.j(xyzTemp[2]);
	    double kTemp = grid.k(xyzTemp[3]);
	    size_t i0 = std::floor(iTemp);  size_t i1 = i0 + 1;
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

	    // Interpolate intensity from neighbouring 8 lattice points to temporary point:
	    Inew[index] = IntegerPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
	    			  intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
	    			  intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1);
    }
	std::swap(I,Inew);
	std::swap(q,qNew);
}



void Radiation::Collide()
{
	PROFILE_FUNCTION();
	// TODO: Steife DGL?
	
	PARALLEL_FOR(3)
	for(size_t k=HALO; k<metric.grid.nz-HALO; k++)
	for(size_t j=HALO; j<metric.grid.ny-HALO; j++)
	for(size_t i=HALO; i<metric.grid.nx-HALO; i++)
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
		
		for(size_t d = 0; d < stencil.nDir; d++)
		{
			size_t index = Index(ijk,d);
			double A = W * (1.0 - Tensor3::Dot(q[ijk] * stencil.Ct3(d), u));

			double Gamma = (eta[ijk] + kappa0[ijk]*fluidE) / (A*A*A) - A*I[index] * (kappaA[ijk] + kappa0[ijk]);
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
		Coord pixel = camera.xyz(ij);   // xyz coord of pixel in world space.
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
	}
}



void Radiation::WriteIntensitiesToCsv(double time, const int frameNumber, std::string directory, std::string name)
{
	PROFILE_FUNCTION();
    CreateDirectory(directory);

    name = name + FrameNumber(frameNumber) + ".csv";
    std::ofstream fileOut(directory + "/" + name);

    fileOut << "#nx=" << grid.nx << "\n";
    fileOut << "#ny=" << grid.ny << "\n";
    fileOut << "#nz=" << grid.nz << "\n";
    fileOut << "#startx=" << grid.startx << "\n";
    fileOut << "#starty=" << grid.starty << "\n";
    fileOut << "#startz=" << grid.startz << "\n";
    fileOut << "#endx=" << grid.endx << "\n";
    fileOut << "#endy=" << grid.endy << "\n";
    fileOut << "#endz=" << grid.endz << "\n";
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
				Tensor3 dir = q[ijk] * stencil.Ct3(d);
                double value = I[index];

				Coord pos = xyz;
				pos[1] += dir[1] * grid.dt * 0.5;// * (0.1 + 0.9 * value);
				pos[2] += dir[2] * grid.dt * 0.5;// * (0.1 + 0.9 * value);
				pos[3] += dir[3] * grid.dt * 0.5;// * (0.1 + 0.9 * value);

				fileOut << pos[1] << ", " << pos[2] << ", " << pos[3] << ", " << value << "\n";
			}
		}
	}
	
    fileOut.close();
}



void Radiation::RunSimulation(Config config)
{
	// -------------------- Initialization --------------------
	int timeSteps = ceil(config.simTime / grid.dt);
	config.simTime = timeSteps * grid.dt;
	Log logger(config.name, config.simTime, stencil, streamingStencil, metric);

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
		std::cout << " nDir         = " << stencil.nDir << "\n";
		std::cout << " nSphHarm     = " << streamingStencil.nDir << "\n";
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
			// Collide();

			if((config.writeData || config.useCamera) && (n % config.writeFrequency) == 0)
				ComputeMomentsLF();
			if(n % config.writeFrequency == 0)
			{
				if(config.writeData)
				{
					grid.WriteFrametoCsv(n*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, n, logger.directoryPath + "/Moments");
					// WriteIntensitiesToCsv(n*grid.dt, n, logger.directoryPath + "/Intensities", "I");
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
			switch(streamingType)
			{
				case(StreamingType::FlatFixed):         StreamFlatFixed();							    break;
				case(StreamingType::FlatAdaptive):      UpdateQuaternions();    StreamFlatAdaptive();	break;
				case(StreamingType::CurvedFixed):       StreamCurvedFixed();	 					    break;
				case(StreamingType::CurvedAdaptive):    UpdateQuaternions();    StreamCurvedAdaptive();	break;
			}

			if(config.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if(config.writeData || config.useCamera)
        {
			ComputeMomentsIF();
			ComputeMomentsLF();
        }
		if(config.writeData)
        {
			grid.WriteFrametoCsv(timeSteps*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, timeSteps, logger.directoryPath + "/Moments");
		    // WriteIntensitiesToCsv(timeSteps*grid.dt, timeSteps, logger.directoryPath + "/Intensities", "I");
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