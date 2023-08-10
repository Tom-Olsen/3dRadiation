#ifndef __INCLUDE_GUARD_Logger_hh__
#define __INCLUDE_GUARD_Logger_hh__
#include <math.h>       // Basic math.
#include "Spacetimes.h" // Metric data.
#include "Stencil.h"    // Velocity stencil.

class Logger
{
public:
    // General simulation parameters:
    Stencil &intensityStencil;
    Stencil &streamingStencil;
    Metric &metric;
    int timeSteps;
    double simTime;
    double interpolationAccuracy;
    double sigmaMax;
    double fluxMax;

    // Data management:
    std::string name;
    std::string date;
    std::string directoryPath;
    std::vector<std::string> timeNames;
    std::vector<double> timeMeasurements;

    Logger(Stencil &intensityStencil, Stencil &streamingStencil, Metric &metric)
        : intensityStencil(intensityStencil), streamingStencil(streamingStencil), metric(metric){};

    void SetValues(std::string name, double simTime, double interpolationAccuracy, double sigmaMax, double fluxMax)
    {
        this->name = name;
        this->interpolationAccuracy = interpolationAccuracy;
        this->sigmaMax = sigmaMax;
        this->fluxMax = fluxMax;

        // Derived from simulation parameters:
        timeSteps = ceil(simTime / metric.grid.dt);
        this->simTime = simTime = timeSteps * metric.grid.dt;

        // Creation time:
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%H.%M.%S - %d.%m.%Y");
        date = oss.str();

        // File system overhead:
        directoryPath = OUTPUTDIR + name;
        if (!CreateDirectory(directoryPath))
            double a = 0;
    }

    ~Logger()
    {
        std::ofstream file(directoryPath + "/Log.txt");

        file << "Creation Date: " << date << std::endl
             << std::endl;

        file << "Spacetime: " << metric.Name() << std::endl;
        file << "Black Hole  Mass,   M = " << metric.m << std::endl;
        file << "Black Hole  Spin,   a = " << metric.a << std::endl
             << std::endl;

        file << "Grid Structure:" << std::endl;
        file << "Simulation Time = " << simTime << std::endl;
        file << "(startx,starty,startz) = (" << metric.grid.startx << "," << metric.grid.starty << "," << metric.grid.startz << ")" << std::endl;
        file << "(endx  ,  endy,  endz) = (" << metric.grid.endx << "," << metric.grid.endy << "," << metric.grid.endz << ")" << std::endl;
        file << "(Nx    ,    Ny,    Nz) = (" << metric.grid.nx << "," << metric.grid.ny << "," << metric.grid.nz << ")" << std::endl;
        file << "Nxyz = " << metric.grid.nxyz << std::endl;
        file << "Nt   = " << timeSteps << std::endl;
        file << "dx   = " << metric.grid.dx << std::endl;
        file << "dy   = " << metric.grid.dy << std::endl;
        file << "dz   = " << metric.grid.dz << std::endl;
        file << "dt   = " << metric.grid.dt << std::endl;
        file << "cfl        = " << metric.grid.GetCFL() << std::endl;
        file << "interp acc = " << 100 * interpolationAccuracy << "%" << std::endl;
        file << "sigma max  = " << sigmaMax << std::endl;
        file << "flux max   = " << fluxMax << std::endl
             << std::endl;

        file << "Stencil Properties:" << std::endl;
        file << "Intensity Stencil: nDir          = " << intensityStencil.nDir << std::endl;
        file << "Intensity Stencil: nOrder        = " << intensityStencil.nOrder << std::endl;
        file << "Intensity Stencil: nCoefficients = " << intensityStencil.nCoefficients << std::endl;
        file << "Streaming Stencil: nDir          = " << streamingStencil.nDir << std::endl;
        file << "Streaming Stencil: nOrder        = " << streamingStencil.nOrder << std::endl;
        file << "Streaming Stencil: nCoefficients = " << streamingStencil.nCoefficients << std::endl
             << std::endl;

        file << "Time measurements:" << std::endl;
        for (int i = 0; i < timeMeasurements.size(); i++)
            file << timeNames[i] << ": " << timeMeasurements[i] << "s" << std::endl;
    }

    void AddTimeMeasurement(std::string timeName, double timeMeasurement)
    {
        timeNames.push_back(timeName);
        timeMeasurements.push_back(timeMeasurement);
    }
};

#endif //__INCLUDE_GUARD_Logger_hh__