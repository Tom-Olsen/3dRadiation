#ifndef __INCLUDE_GUARD_Camera_h__
#define __INCLUDE_GUARD_Camera_h__
#include "TensorTypes.hh"



class Camera
{
public:
    double* image;
    int resX;
    int resY;
    int pixelCount;
    int width;
    int height;

    Coord position;
    glm::vec3 eulerAngles;
    Tensor3 lookDirection;
    Tensor3 orthogonalPassThrough;

    Camera(int resX, int resY, int width, int height, Coord position, glm::vec3 eulerAngles);
    Camera();
    
    int Index(int i, int j);
    Coord xyz(int i, int j);
    Coord xyz(int ij);

    void WriteImagetoCsv
    (float time, const int frameNumber, std::string directory, std::string name="");
};



#endif //__INCLUDE_GUARD_Camera_h__