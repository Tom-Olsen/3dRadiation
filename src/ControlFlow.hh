#ifndef __INCLUDE_GUARD_ControlFlow_hh__
#define __INCLUDE_GUARD_ControlFlow_hh__



// Frames:
class IF {};
class LF {};



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