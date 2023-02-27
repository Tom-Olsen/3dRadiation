#include "Camera.h"



Camera::Camera(int resX, int resY, int width, int height, Coord position, glm::vec3 eulerAngles) :
resX(resX), resY(resY), width(width), height(height), position(position), eulerAngles(eulerAngles)
{
    pixelCount = resX * resY;
    image = new double[pixelCount];

    glm::quat q(eulerAngles);
    lookDirection = q * Tensor3(0,0,1);
    orthogonalPassThrough[1] = -lookDirection[1];
    orthogonalPassThrough[2] = -lookDirection[2];
    orthogonalPassThrough[3] = -lookDirection[3];
}
Camera::Camera()
{
    resX = resY = 10;
    width = height = 1;
    position = Coord(0.0);
    eulerAngles = glm::vec3(0,0,0);

    pixelCount = resX * resY;
    image = new double[pixelCount];

    glm::quat q(eulerAngles);
    lookDirection = q * Tensor3(0,0,1);
    orthogonalPassThrough[1] = -lookDirection[1];
    orthogonalPassThrough[2] = -lookDirection[2];
    orthogonalPassThrough[3] = -lookDirection[3];
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
Coord Camera::xyz(int ij)
{
    int i = ij % resX;
    int j = ij / resX;
    return xyz(i,j);
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
        fileOut << x[1] << "," << x[2] << "," << x[3] << "," << image[ij] << "\n";
    }

    // Look Direction:
    int nPoints = 20;
    for(int i=0; i<nPoints; i++)
    {
        double dist = sqrt(width * height) / 5.0 * i / (double)nPoints;
        fileOut << position[1] + dist * lookDirection[1] << "," << position[2] + dist * lookDirection[2] << "," << position[3] + dist * lookDirection[3] << "," << 0 << "\n";
    }


    fileOut.close();
}