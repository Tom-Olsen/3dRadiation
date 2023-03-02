#ifndef __INCLUDE_GUARD_Spacetimes_h__
#define __INCLUDE_GUARD_Spacetimes_h__
#include <string>
#include "ControlFlow.hh"   // Template arguments and profiling macros.
#include "TensorTypes.hh"   // General relativity tensors.
#include "Grid.h"           // Numerical Grid and mapping to physical domain.
#include "Metric.h"         // Metric parent class.



class Minkowski : public Metric
{
public:
    Minkowski(Grid& grid_, double m_=0, double a_=0);
    Minkowski(const Minkowski& minkowski) = delete;

    bool InsideBH(const Coord& xyz) override;
    std::string Name() override;
private:
    Tensor4x4 MetricFunction(const Coord& xyz) override;
};



class SchwarzSchild : public Metric
{
public:
    SchwarzSchild(Grid& grid_, double m_=1, double a_=0);
    SchwarzSchild(const SchwarzSchild& schwarzSchild) = delete;

    bool InsideBH(const Coord& xyz) override;
    std::string Name() override;
private:
    Tensor4x4 MetricFunction(const Coord& xyz) override;
};



class KerrSchild : public Metric
{
public:
    KerrSchild(Grid& grid_, double m_=1, double a_=0);
    KerrSchild(const KerrSchild& kerrSchild) = delete;

    bool InsideBH(const Coord& xyz) override;
    std::string Name() override;
private:
    Tensor4x4 MetricFunction(const Coord& xyz) override;
};
#endif //__INCLUDE_GUARD_Spacetimes_h__