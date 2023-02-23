#include "Camera.h"



Camera::Camera(int resX, int resY, int width, int height, Coord position, glm::vec3 eulerAngles) :
resX(resX), resY(resY), width(width), height(height), position(position), eulerAngles(eulerAngles)
{
    image = new double[resX * resY];
}

Camera::Camera(const Camera& camera) :
resX(camera.resX), resY(camera.resY), width(camera.width), height(camera.height),
position(camera.position), eulerAngles(camera.eulerAngles)
{
    image = new double[resX * resY];
}


int Camera::Index(int i, int j)
{ return i + j * resX; }

Coord Camera::xyz(int i, int j)
{
    // xLocal € [ -width/2,  width/2]
    // yLocal € [-height/2, height/2]
    double xLocal = width  * (i / (resX - 1.0) - 0.5);
    double yLocal = height * (j / (resY - 1.0) - 0.5);
    double zLocal = 0;
    Coord xyzLocal(xLocal,yLocal,zLocal);

    glm::quat q(eulerAngles);

    Coord xyz = q * xyzLocal + position;
    return xyz;
}

void Camera::WriteImagetoCsv
(float time, const int frameNumber, std::string directory, std::string name)
{
    CreateDirectory(directory);
    
    name = (name == "") ? "image" :  name;
    std::ofstream fileOut(directory + "/" + name + FrameNumber(frameNumber) + "_" + std::to_string(resX) + "x" + std::to_string(resY) + "y" + ".csv");
    fileOut << "#x,y,z,value\n";

    for(int j=0; j<resY; j++)
    for(int i=0; i<resX; i++)
    {
        int ij = Index(i,j);
        Coord x = xyz(i,j);
        fileOut << x[1]   << "," << x[2]   << "," << x[3] << "," << image[ij] << "\n";
    }

    fileOut.close();
}