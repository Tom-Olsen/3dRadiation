#ifndef __INCLUDE_GUARD_ControlFlow_hh__
#define __INCLUDE_GUARD_ControlFlow_hh__

// Frames:
class IF
{
};
class LF
{
};

// If defined std::sin etc are used. Otherwise the Mysin etc. optimisations are used.
// #define USE_STD_MATH
// Order of polynomial approximation for MySin, MyCos, MyAtan, MyAtan2:
#define APPROXIMATION_ORDER 9

// Output Directory:
#define OUTPUTDIR (std::string) "../output/"

// omp parallel for macro:
#define STRINGIFY(X) #X
#define PRAGMA(X) _Pragma(STRINGIFY(X))
// #define PARALLEL_FOR(N)
#define PARALLEL_FOR(N) PRAGMA(omp parallel for collapse(N) schedule(dynamic,50))

// Time measurements:
#define PROFILING 1
#if PROFILING
#define PROFILE_SCOPE(name) Profiler::Timer timer##__LINE__(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__PRETTY_FUNCTION__)
#else
#define PROFILE_SCOPE(name)
#define PROFILE_FUNCTION()
#endif

#endif //__INCLUDE_GUARD_ControlFlow_hh__