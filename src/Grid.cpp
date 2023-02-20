#include "Grid.h"



// Constructors:
Grid::Grid(int nx, int ny, int nz, Coord start, Coord end) :
nx(nx), ny(ny), nz(nz),
startx(start[1]), starty(start[2]), startz(start[3]), 
endx(end[1]), endy(end[2]), endz(end[3])
{
    if(end[1] < start[1])
        exit_on_error("Grid x€[a,b] has b<a!");
    if(end[2] < start[2])
        exit_on_error("Grid y€[a,b] has b<a!");
    if(end[3] < start[3])
        exit_on_error("Grid z€[a,b] has b<a!");
    if(nx==1 or ny==1 or nz==1)
        exit_on_error("Grid must have at least 2 LP in each Dimension.");

    dx = (endx - startx) / (nx - 1.0);
    dy = (endy - starty) / (ny - 1.0);
    dz = (endz - startz) / (nz - 1.0);
    dt = m_cfl * std::min({dx, dy, dz});
    nxy  = nx * ny;
    nxyz = nx * ny * nz;
}
Grid::Grid(const Grid& grid) :
nx(grid.nx), ny(grid.ny), nz(grid.nz),
startx(grid.startx), starty(grid.starty), startz(grid.startz), 
endx(grid.endx), endy(grid.endy), endz(grid.endz)
{
    dx = (endx - startx) / (nx - 1);
    dy = (endy - starty) / (ny - 1);
    dz = (endz - startz) / (nz - 1);
    dt = m_cfl * std::min({dx, dy, dz});
    nxy  = nx * ny;
    nxyz = nx * ny * nz;
}

// Setters/Getters:
void Grid::SetCFL(double cfl)
{
    m_cfl = cfl;
    dt = m_cfl * std::min({dx, dy, dz});
}
double Grid::GetCFL()
{ return m_cfl; }

// Grid Access Tools:
int Grid::Index(int i, int j, int k)
{ return i + j * nx + k * nxy; }
int Grid::Index(int i, int j, int k, int d)
{ return i + j * nx + k * nxy + d * nxyz; }
Coord Grid::xyz(int i, int j, int k)
{ return Coord(startx + i*dx, starty + j*dy, startz + k*dz); }
Coord Grid::xyz(double i, double j, double k)
{ return Coord(startx + i*dx, starty + j*dy, startz + k*dz); }
double Grid::i(double x)
{ return (x - startx) / dx; }
double Grid::j(double y)
{ return (y - starty) / dy; }
double Grid::k(double z)
{ return (z - startz) / dz; }
Coord Grid::ijk(const Coord& xyz)
{ return Coord(i(xyz[1]), j(xyz[2]), k(xyz[3])); }

// Domain Checks:
bool Grid::OutsideDomain(const Coord& xyz)
{
    return
    xyz[1]>endx or xyz[1]<startx or
    xyz[2]>endy or xyz[2]<starty or
    xyz[3]>endz or xyz[3]<startz;
}
bool Grid::OutsideDomain(double i, double j, double k)
{
    return
    i > nx - 1 or i < 0 or
    j > ny - 1 or j < 0 or
    k > nz - 1 or k < 0;
}

// Write Data to file:
void Grid::WriteFrametoJson
(float time, double* r, double* g, double* b, double* a,
 const int frameNumber, std::string directory, std::string name)
{
    PROFILE_FUNCTION();
    // main body:
    Json::Value jsonData;

    // primitive types:
    jsonData["meshType"] = 0;
    jsonData["time"] = time;

    // arrays:
    Json::Value colors(Json::arrayValue);
    for(int ijk=0; ijk<nxyz; ijk++)
    {
        Json::Value Color;
        Color["r"] = r[ijk];
        Color["g"] = g[ijk];
        Color["b"] = b[ijk];
        Color["a"] = a[ijk];
        colors.append(Color);
    }
    jsonData["colors"] = colors;

    // structs:
    Json::Value start;
    start["x"] = startx;
    start["y"] = starty;
    start["z"] = startz;
    jsonData["start"] = start;
    
    Json::Value end;
    end["x"] = endx;
    end["y"] = endy;
    end["z"] = endz;
    jsonData["end"] = end;
    
    Json::Value resolution;
    resolution["x"] = nx;
    resolution["y"] = ny;
    resolution["z"] = nz;
    jsonData["resolution"] = resolution;
    
    // write json to file:
    // TODO: check if directory exists and dont create it if not necessary
    CreateDirectory(directory);
    name = (name == "") ? "data" :  name;

    std::ofstream fileOut(directory + "/" + name + FrameNumber(frameNumber) + ".json");
    fileOut << jsonData;
    fileOut.close();
}
void Grid::WriteFrametoCsv
(float time, double* r, double* g, double* b, double* a,
 const int frameNumber, std::string directory, std::string name)
{
    CreateDirectory(directory);
    
    name = (name == "") ? "data" :  name;
    std::ofstream fileOut(directory + "/" + name + FrameNumber(frameNumber) + "_" + std::to_string(nx) + "x" + std::to_string(ny) + "y" + std::to_string(nz) + "z" +  ".csv");
    fileOut << "#x,y,z,r,g,b,a\n";

    for(int k=0; k<nz; k++)
    for(int j=0; j<ny; j++)
    for(int i=0; i<nx; i++)
    {
        int ijk = Index(i,j,k);
        Coord x = xyz(i,j,k);
        fileOut << x[1]   << "," << x[2]   << "," << x[3] << ",";
        fileOut << r[ijk] << "," << g[ijk] << "," << b[ijk] << "," << a[ijk] << "\n";
    }

    fileOut.close();
}