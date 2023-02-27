#ifndef __INCLUDE_GUARD_Utility_hh__
#define __INCLUDE_GUARD_Utility_hh__
#include <iomanip>              // std::setprecision(), std::put:time
#include <iostream>             // Output to terminal.
#include <fstream>              // File input/output.
#include <cmath>                // signbit.
#include <algorithm>            // clamp.
#include <sstream>              // stringstream.
#include <filesystem>           // Folder/file management.
#include <vector>               // basic c++ vector
#include "ControlFlow.hh"       // Template arguments and profiling macros.
#include "eigen/Eigen/Dense"    // Eigen library for solving linear systems.



// Felix:
template <class T>
class AlignedAllocator
{
public:
    using value_type = T;

    // Constructor:
    T * allocate(size_t size)
    { return static_cast<T *>(std::aligned_alloc(64,size * sizeof(T))); }

    // Descturctor:
    void deallocate(T * p_t, size_t)
    { std::free(p_t); }

    // Initialize each element:
    template<class U, class... Args>
    void construct(U* p, Args&&... args)
    { *p = 0; }

    // Needed for std::swap(..., ...):
    bool operator==(const AlignedAllocator &)
    { return true; }
};
using RealBuffer = std::vector<double,AlignedAllocator<double>>;



// Fast integer exponentiation.
template<int N>
inline double MyPow(double a)
{ return a * MyPow<N-1>(a); }
template<>
inline double MyPow<0>(double a)
{ return 1; }

// Get sign of T
template <typename T>
inline int sgn(T val)
{ return (T(0) < val) - (val < T(0)); }



// Marker for debugging
inline void Marker(std::string name="", bool newline=true)
{
    static int i = 0;
    std::cout << name << ": " << i << std::endl;
    if(newline)
        std::cout << std::endl;
    i++;
}



// Number of digits of given number.
inline int numDigits(int number)
{
    int digits = 0;
    while (number)
	{
        number /= 10;
        digits++;
    }
    return digits;
}



// Converts frame number to correct format, e.g. 1 -> 0001.
inline std::string FrameNumber(unsigned int f)
{
	std::string frameNumber;
	int maxDigits = 4;
	int fDigits   = numDigits(f);

	for(int i=0; i<(maxDigits-fDigits); i++)
		frameNumber += "0";
    if(f!=0)
    	frameNumber += std::to_string(f);

	return frameNumber;
}



inline std::string Format(const double n, const int precision=6)
{
    std::string output;

    // Leading space for positive numbers:
    if(n == 0 and std::signbit(n) == 0){output = " ";}
    else if(n == 0 and std::signbit(n) == 1){output = "";}
	else if(n >= 0){output = " ";}
	else   {output = "";}

    // Leading space if 2 digit number:
    //if(abs(n)<100){output += " ";}
    // Leading space if 1 digit number:
    //if(abs(n)<10){output += " ";}

    // Number of decimal digits:
    std::ostringstream ss;
    ss.precision(precision);
    ss << std::fixed << n; 
    output += ss.str();

    return output;    
}



// Print formatted double.
inline void PrintDouble(const double d, std::string name, bool newline=false, const int precision=6)
{
    std::cout << name << " = "<< Format(d,precision) << "\n";
    if(newline) std::cout << "\n";
}



// Exit programm with Error Message.
inline void exit_on_error(const char* const msg="")
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(errno);
}



// Print formatted Matrix.
inline void PrintMat(const Eigen::MatrixXd& A, int precision=6)
{
    for(int i=0; i<A.cols(); i++)
    {
        std::cout << "(" << Format(A(i,0),precision);
        for(int j=1; j<A.rows(); j++)
            std::cout << ", " << Format(A(i,j),precision);
        std::cout << ")" << std::endl;
    }
}



// File management:
// Creates directory in current working directory.
// Nested directories possible, e.g: path="data/output/xValues"
inline bool CreateDirectory(std::string path)
{
	namespace fs = std::filesystem;
    fs::path currentPath = fs::current_path();
	std::string directoryPath = currentPath/path;
    return fs::create_directories(directoryPath);
}
#endif //__INCLUDE_GUARD_Utility_hh__