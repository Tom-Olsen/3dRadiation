#ifndef __INCLUDE_GUARD_Camera_h__
#define __INCLUDE_GUARD_Camera_h__
#include "TensorTypes.hh"



class Camera
{
public:
    double* image;
    int resX;
    int resY;
    int width;
    int height;

    Coord position;
    glm::vec3 eulerAngles;

    Camera(int resX, int resY, int width, int height, Coord position, glm::vec3 eulerAngles);
    Camera() = delete;
    Camera(const Camera& Camera);
    
    int Index(int i, int j);
    Coord xyz(int i, int j);

    void WriteImagetoCsv
    (float time, const int frameNumber, std::string directory, std::string name="");
};



#endif //__INCLUDE_GUARD_Camera_h__