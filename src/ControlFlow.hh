#ifndef __INCLUDE_GUARD_ControlFlow_hh__
#define __INCLUDE_GUARD_ControlFlow_hh__



// Frames:
class IF {};
class LF {};



// Datastructure:
// #define ijkd0d1
#define d0d1ijk // is faster on CPU



// Order of polynomial approximation for MySin, MyCos, MyAtan, MyAtan2:
#define APPROXIMATION_ORDER 9



// IntensityTypes:
class Bulk {};
class North {};
class South {};



// StaticOrDynamic:
class Static {};
class Dynamic {};



// Output Directory:
// #define OUTPUTDIR "output/"
#define OUTPUTDIR "/mnt/ceph/tolsen/output/"



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