#ifndef __INCLUDE_GUARD_Camera_h__
#define __INCLUDE_GUARD_Camera_h__
#include "Utility.hh"                   // Utility functions.
#include "Profiler.hh"                  // Time measurement profiler.
#include "DataTypes.hh" 			    // General relativity tensors.
#include "glm/glm/gtc/quaternion.hpp"   // Quaternions.



class Camera
{
public:
    RealBuffer image;
    size_t resX;
    size_t resY;
    size_t pixelCount;
    size_t width;
    size_t height;

    Coord position;
    glm::vec3 eulerAngles;
    Tensor3 lookDirection;

    Camera(size_t resX, size_t resY, size_t width, size_t height, Coord position, glm::vec3 eulerAngles);
    Camera();
    
    size_t Index(size_t i, size_t j);
    Coord xyzLocal(size_t i, size_t j);
    Coord xyzLocal(size_t ij);
    Coord xyzWorld(size_t i, size_t j);
    Coord xyzWorld(size_t ij);


    void WriteImagetoCsv(float time, std::string directory, std::string name = "");
};



#endif //__INCLUDE_GUARD_Camera_h__