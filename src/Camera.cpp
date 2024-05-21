#include "Camera.h"



Camera::Camera(size_t resX, size_t resY, size_t width, size_t height, Coord position, glm::vec3 eulerAngles) :
resX(resX), resY(resY), width(width), height(height), position(position), eulerAngles(eulerAngles)
{
    pixelCount = resX * resY;
    image.resize(pixelCount);

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
    image.resize(pixelCount);

    glm::quat q(eulerAngles);
    lookDirection = q * Tensor3(0, 0, 1);
}



size_t Camera::Index(size_t i, size_t j)
{ return i + j * resX; }

Coord Camera::xyzLocal(size_t i, size_t j)
{
    double x = width  * (i / (resX - 1.0) - 0.5); // xLocal € [ -width/2,  width/2]
    double y = height * (j / (resY - 1.0) - 0.5); // yLocal € [-height/2, height/2]
    double z = 0;
    return Coord(x, y, z);
}
Coord Camera::xyzLocal(size_t ij)
{
    size_t i = ij % resX;
    size_t j = ij / resX;
    return xyzLocal(i, j);
}
Coord Camera::xyzWorld(size_t i, size_t j)
{
    Coord xyz = glm::quat(eulerAngles) * xyzLocal(i, j) + position;
    return xyz;
}
Coord Camera::xyzWorld(size_t ij)
{
    size_t i = ij % resX;
    size_t j = ij / resX;
    return xyzWorld(i, j);
}

void Camera::WriteImagetoCsv(float time, std::string directory, std::string name)
{
	PROFILE_FUNCTION();
    CreateDirectory(directory);
    
    name = (name == "") ? "image" :  name;
    std::ofstream fileOut(directory + name + Format(time, 3) + "t " + std::to_string(resX) + "nx " + std::to_string(resY) + "ny" + ".csv");
    
    fileOut << "#nx=" << resX << "\n";
    fileOut << "#ny=" << resY << "\n";
    fileOut << "#startx=" << position[0] - width  / 2.0 << "\n";
    fileOut << "#starty=" << position[1] - height / 2.0 << "\n";
    fileOut << "#endx=" << position[0] + width  / 2.0 << "\n";
    fileOut << "#endy=" << position[1] + height / 2.0 << "\n";
    fileOut << "#eulerx=" << eulerAngles[0] << "\n";
    fileOut << "#eulery=" << eulerAngles[1] << "\n";
    fileOut << "#eulerz=" << eulerAngles[2] << "\n";
    fileOut << "#t=" << time << "\n";
    fileOut << "#x,y,r,g,b,a\n";

    for(size_t j=0; j<resY; j++)
        for(size_t i=0; i<resX; i++)
        {
            size_t ij = Index(i ,j);
            Coord x = xyzLocal(i, j);
            fileOut << x[1] << "," << x[2] << ",";
            fileOut << image[ij] << "," << 0 << "," << 0 << "," << 0 << "\n";
        }

    fileOut.close();
}