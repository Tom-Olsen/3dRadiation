#include "Radiation2.h"



Radiation2::Radiation2(Metric& metric, Stencil& intensityStencil, Stencil& streamingStencil, Camera& camera):
grid(metric.grid), metric(metric), intensityStencil(intensityStencil), streamingStencil(streamingStencil), camera(camera)
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
	coefficientsI.resize(grid.nxyz * intensityStencil.nCoefficients);
	coefficientsInew.resize(grid.nxyz * intensityStencil.nCoefficients);
	coefficientsS.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsX.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsY.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsZ.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsCx.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsCy.resize(grid.nxyz * streamingStencil.nCoefficients);
	coefficientsCz.resize(grid.nxyz * streamingStencil.nCoefficients);
}
Radiation2::~Radiation2()
{
	delete[] isInitialGridPoint;
}




size_t Radiation2::IntensityHarmonicIndex(size_t ijk, size_t d)
{
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return d + ijk * intensityStencil.nCoefficients;
	#endif
}
size_t Radiation2::IntensityHarmonicIndex(size_t i, size_t j, size_t k, size_t d)
{
	size_t ijk = grid.Index(i,j,k);
	#ifdef ijkd
		return ijk + d * grid.nxyz;
	#endif
	#ifdef dijk
		return d + ijk * intensityStencil.nCoefficients;
	#endif
}
size_t Radiation2::StreamingHarmonicIndex(size_t ijk, size_t d)
{
	return d + ijk * streamingStencil.nCoefficients;
}



void Radiation2::NormalizeInitialDirections()
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



void Radiation2::LoadInitialData()
{
	PROFILE_FUNCTION();
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
                double I[intensityStencil.nDir];
				for(size_t d=0; d<intensityStencil.nDir; d++)
					I[d] = initialE[ijk];
                SphericalHarmonicsXyz::GetCoefficients(intensityStencil, I, &coefficientsI[IntensityHarmonicIndex(ijk,0)]);
			}
			else
			{// Kent intensity distribution: https://en.wikipedia.org/wiki/Kent_distribution
			 	Tensor3 n = Tensor3(initialNx[ijk],initialNy[ijk],initialNz[ijk]);
                double I[intensityStencil.nDir];
				for(size_t d=0; d<intensityStencil.nDir; d++)
				{
					Tensor3 p = intensityStencil.Ct3(d);
					I[d] = initialE[ijk] * exp(sigma * Tensor3::Dot(n, p));
				}
                SphericalHarmonicsXyz::GetCoefficients(intensityStencil, I, &coefficientsI[IntensityHarmonicIndex(ijk,0)]);
			}
		}
	}
}



void Radiation2::NormalizeInitialIntensities()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ijk=0; ijk<grid.nxyz; ijk++)
		if(isInitialGridPoint[ijk])
			initialE[ijk] *= initialE[ijk] / E_LF[ijk];
}



void Radiation2::UpdateStreamingCoefficients()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
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
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataS , &coefficientsS[StreamingHarmonicIndex(ijk,0)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataX , &coefficientsX[StreamingHarmonicIndex(ijk,0)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataY , &coefficientsY[StreamingHarmonicIndex(ijk,0)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataZ , &coefficientsZ[StreamingHarmonicIndex(ijk,0)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[StreamingHarmonicIndex(ijk,0)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[StreamingHarmonicIndex(ijk,0)]);
		SphericalHarmonicsXyz::GetCoefficients(streamingStencil, dataCz, &coefficientsCz[StreamingHarmonicIndex(ijk,0)]);
	}
}



void Radiation2::ComputeMomentsIF()
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
		for(size_t d=0; d<intensityStencil.nDir; d++)
		{
			Tensor3 dir = intensityStencil.Ct3(d);
            double I_d = intensityStencil.W(d) * IntensityAt(ijk,dir);
			E[ijk]   += I_d;
			Fx[ijk]  += I_d * dir[1];
			Fy[ijk]  += I_d * dir[2];
			Fz[ijk]  += I_d * dir[3];
			Pxx[ijk] += I_d * dir[1] * dir[1];
			Pxy[ijk] += I_d * dir[1] * dir[2];
			Pxz[ijk] += I_d * dir[1] * dir[3];
			Pyy[ijk] += I_d * dir[2] * dir[2];
			Pyz[ijk] += I_d * dir[2] * dir[3];
			Pzz[ijk] += I_d * dir[3] * dir[3];
		}
        // E[ijk]   *= fourPiInv;
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
void Radiation2::ComputeMomentsLF()
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



Coord Radiation2::GetTempCoordinate(size_t ijk, Tensor3 direction)
{
	Coord xyzTemp;
	xyzTemp[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsX[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
	xyzTemp[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsY[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
	xyzTemp[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsZ[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
	return xyzTemp;
}
Tensor3 Radiation2::GetTemp3VelocityIF(size_t ijk, Tensor3 direction)
{
	Tensor3 vTempIF;
	vTempIF[1] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCx[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
	vTempIF[2] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCy[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
	vTempIF[3] = SphericalHarmonicsXyz::GetValue(direction, &coefficientsCz[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
	return vTempIF;
}
double Radiation2::GetFrequencyShift(size_t ijk, Tensor3 direction)
{
	return SphericalHarmonicsXyz::GetValue(direction, &coefficientsS[StreamingHarmonicIndex(ijk,0)], streamingStencil.nCoefficients);
}
double Radiation2::IntensityAt(size_t ijk, const Tensor3& vTempIF)
{
	return SphericalHarmonicsXyz::GetValue(vTempIF, &coefficientsI[IntensityHarmonicIndex(ijk,0)], intensityStencil.nCoefficients);
}



void Radiation2::Stream()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	for(size_t k=2; k<grid.nz-2; k++)
	for(size_t j=2; j<grid.ny-2; j++)
	for(size_t i=2; i<grid.nx-2; i++)
	{
		int ijk = grid.Index(i,j,k);
		double alpha = metric.GetAlpha(ijk);
		double I[intensityStencil.nDir];
		for(size_t d=0; d<intensityStencil.nDir; d++)
		{
			// Skip LPs which are inside BH:
			if(metric.InsideBH(grid.xyz(i,j,k)))
			{
				I[d] = 0;
				continue;
			}

			// Get quantities at emission point:
			Tensor3 direction = intensityStencil.Ct3(d);
			double s = GetFrequencyShift(ijk, direction);
			Coord xyzTemp = GetTempCoordinate(ijk, direction);
			Tensor3 vTempIF = GetTemp3VelocityIF(ijk, direction);
			
			if(metric.InsideBH(xyzTemp))
			{
				I[d] = 0;
				continue;
			}

			// Get 8 nearest Grid Points:
			double iTemp = grid.i(xyzTemp[1]);
			double jTemp = grid.j(xyzTemp[2]);
			double kTemp = grid.k(xyzTemp[3]);
			size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
			size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
			size_t k0 = std::floor(kTemp);	size_t k1 = k0 + 1;
			
			// Bilinear intensity interpolation in velocity space:
			double intensityAt_i0j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k0))) * IntensityAt(grid.Index(i0,j0,k0),vTempIF);
			double intensityAt_i0j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0,k1))) * IntensityAt(grid.Index(i0,j0,k1),vTempIF);
			double intensityAt_i0j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k0))) * IntensityAt(grid.Index(i0,j1,k0),vTempIF);
			double intensityAt_i0j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1,k1))) * IntensityAt(grid.Index(i0,j1,k1),vTempIF);
			double intensityAt_i1j0k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k0))) * IntensityAt(grid.Index(i1,j0,k0),vTempIF);
			double intensityAt_i1j0k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0,k1))) * IntensityAt(grid.Index(i1,j0,k1),vTempIF);
			double intensityAt_i1j1k0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k0))) * IntensityAt(grid.Index(i1,j1,k0),vTempIF);
			double intensityAt_i1j1k1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1,k1))) * IntensityAt(grid.Index(i1,j1,k1),vTempIF);

			// Trilinear intensity interpolation in grid space:
			I[d] = std::max(0.0,
			 IntegerPow<4>(s) * TrilinearInterpolation(iTemp-i0, jTemp-j0, kTemp-k0,
			 intensityAt_i0j0k0, intensityAt_i0j0k1, intensityAt_i0j1k0, intensityAt_i0j1k1,
			 intensityAt_i1j0k0, intensityAt_i1j0k1, intensityAt_i1j1k0, intensityAt_i1j1k1));
		}
		SphericalHarmonicsXyz::GetCoefficients(intensityStencil, I, &coefficientsInew[IntensityHarmonicIndex(ijk,0)]);
	}
	std::swap(coefficientsI, coefficientsInew);
}



void Radiation2::RunSimulation(Config config)
{
	// -------------------- Initialization --------------------
	int timeSteps = ceil(config.simTime / grid.dt);
	config.simTime = timeSteps * grid.dt;
	Log2 logger(config.name, config.simTime, intensityStencil, streamingStencil, metric);

	Profiler::Session& session = Profiler::Session::Get();
	session.Start(config.name, "output/" + config.name + "/profileResults.json");
	
	NormalizeInitialDirections();
	LoadInitialData();
	ComputeMomentsIF();
	ComputeMomentsLF();
	NormalizeInitialIntensities();
	LoadInitialData(); // loads normalized data
	UpdateStreamingCoefficients();
	
	// Initial data output:
	if (config.printToTerminal)
	{
		std::cout << " sigma         = " << sigma << "\n";
		std::cout << " nx            = " << grid.nx << "\n";
		std::cout << " ny            = " << grid.ny << "\n";
		std::cout << " nz            = " << grid.nz << "\n";
		std::cout << " nDirStreaming = " << streamingStencil.nDir << "\n";
		std::cout << " nDirIntensity = " << intensityStencil.nDir << "\n";
		std::cout << " simTime       = " << config.simTime << "\n";
		std::cout << " dx            = " << grid.dx << "\n";
		std::cout << " dy            = " << grid.dy << "\n";
		std::cout << " dz            = " << grid.dz << "\n";
		std::cout << " dt            = " << grid.dt << "\n";
		std::cout << " timeSteps     = " << timeSteps << "\n";
		std::cout << " filesToWrite  = " << timeSteps / config.writeFrequency << std::endl;
	}
	// --------------------------------------------------------



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
			if(config.writeData && (n % config.writeFrequency) == 0)
				grid.WriteFrametoCsv(n*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, n, logger.directoryPath + "/Moments", config.name);
			//if(config.useCamera && (n % config.writeFrequency) == 0)
			//{
			//	TakePicture();
			//	camera.WriteImagetoCsv(n*grid.dt, n, logger.directoryPath + "/CameraImages");
			//}

			if (config.updateSphericalHarmonics)
			{ UpdateStreamingCoefficients(); }

			Stream();

			if(config.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if(config.writeData || config.useCamera)
			ComputeMomentsLF();
		if(config.writeData)
			grid.WriteFrametoCsv(timeSteps*grid.dt, E_LF, Fx_LF, Fy_LF, Fz_LF, timeSteps, logger.directoryPath + "/Moments", config.name);
		//if(config.useCamera)
		//{
		//	TakePicture();
		//	camera.WriteImagetoCsv(timeSteps*grid.dt, timeSteps, logger.directoryPath + "/CameraImages");
		//}
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