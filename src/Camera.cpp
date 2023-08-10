#include "Camera.h"



Camera::Camera(int resX, int resY, int width, int height, Coord position, glm::vec3 eulerAngles) :
resX(resX), resY(resY), width(width), height(height), position(position), eulerAngles(eulerAngles)
{
    pixelCount = resX * resY;
    image = new double[pixelCount];

    glm::quat q(eulerAngles);
    lookDirection = q * Tensor3(0,0,1);
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
}



size_t Camera::Index(size_t i, size_t j)
{ return i + j * resX; }

Coord Camera::xyz(size_t i, size_t j)
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
Coord Camera::xyz(size_t ij)
{
    size_t i = ij % resX;
    size_t j = ij / resX;
    return xyz(i,j);
}

void Camera::WriteImagetoCsv
(float time, const int frameNumber, std::string directory, std::string name)
{
	PROFILE_FUNCTION();
    CreateDirectory(directory);
    
    name = (name == "") ? "image" :  name;
    std::ofstream fileOut(directory + name + FrameNumber(frameNumber) + "_" + std::to_string(resX) + "x" + std::to_string(resY) + "y" + ".csv");
    
    fileOut << "#nx=" << resX << "\n";
    fileOut << "#ny=" << resY << "\n";
    fileOut << "#nz=" << 1 << "\n";
    fileOut << "#startx=" << position[0] - width  / 2.0 << "\n";
    fileOut << "#starty=" << position[1] - height / 2.0 << "\n";
    fileOut << "#startz=" << -0.001 << "\n";
    fileOut << "#endx=" << position[0] + width  / 2.0 << "\n";
    fileOut << "#endy=" << position[1] + height / 2.0 << "\n";
    fileOut << "#endz=" << 0.001 << "\n";
    fileOut << "#eulerx=" << eulerAngles[0] << "\n";
    fileOut << "#eulery=" << eulerAngles[1] << "\n";
    fileOut << "#eulerz=" << eulerAngles[2] << "\n";
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