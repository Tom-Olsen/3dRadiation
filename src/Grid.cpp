#include "Grid.h"

// Constructors:
Grid::Grid(size_t nx, size_t ny, size_t nz, Coord start, Coord end) : nx(nx), ny(ny), nz(nz),
                                                                      startx(start[1]), starty(start[2]), startz(start[3]),
                                                                      endx(end[1]), endy(end[2]), endz(end[3])
{
    if (end[1] < start[1])
        ExitOnError("Grid x€[a,b] has b<a!");
    if (end[2] < start[2])
        ExitOnError("Grid y€[a,b] has b<a!");
    if (end[3] < start[3])
        ExitOnError("Grid z€[a,b] has b<a!");
    if (nx == 1 or ny == 1 or nz == 1)
        ExitOnError("Grid must have at least 2 LP in each Dimension.");

    dx = (endx - startx) / (nx - 1.0);
    dy = (endy - starty) / (ny - 1.0);
    dz = (endz - startz) / (nz - 1.0);
    dxInv = 1.0 / dx;
    dyInv = 1.0 / dy;
    dzInv = 1.0 / dz;
    dt = m_cfl * std::min({dx, dy, dz});
    nxy = nx * ny;
    nxyz = nx * ny * nz;
}
Grid::Grid(const Grid &grid) : nx(grid.nx), ny(grid.ny), nz(grid.nz),
                               startx(grid.startx), starty(grid.starty), startz(grid.startz),
                               endx(grid.endx), endy(grid.endy), endz(grid.endz)
{
    dx = (endx - startx) / (nx - 1);
    dy = (endy - starty) / (ny - 1);
    dz = (endz - startz) / (nz - 1);
    dt = m_cfl * std::min({dx, dy, dz});
    nxy = nx * ny;
    nxyz = nx * ny * nz;
}

// Setters/Getters:
void Grid::SetCFL(double cfl)
{
    m_cfl = cfl;
    dt = m_cfl * std::min({dx, dy, dz});
}
double Grid::GetCFL()
{
    return m_cfl;
}

// Grid Access Tools:
size_t Grid::Index(size_t i, size_t j, size_t k)
{
    return i + j * nx + k * nxy;
}
double Grid::x(size_t i)
{
    return startx + i * dx;
}
double Grid::y(size_t j)
{
    return starty + j * dy;
}
double Grid::z(size_t k)
{
    return startz + k * dz;
}
double Grid::x(double i)
{
    return startx + i * dx;
}
double Grid::y(double j)
{
    return starty + j * dy;
}
double Grid::z(double k)
{
    return startz + k * dz;
}
Coord Grid::xyz(size_t i, size_t j, size_t k)
{
    return Coord(x(i), y(j), z(k));
}
Coord Grid::xyz(double i, double j, double k)
{
    return Coord(x(i), y(j), z(k));
}
Coord Grid::xyz(size_t ijk)
{
    return xyz(i(ijk), j(ijk), k(ijk));
}
double Grid::i(double x)
{
    return (x - startx) * dxInv;
}
double Grid::j(double y)
{
    return (y - starty) * dyInv;
}
double Grid::k(double z)
{
    return (z - startz) * dzInv;
}
size_t Grid::i(size_t ijk)
{
    size_t K = k(ijk);
    size_t J = j(ijk);
    return ijk - J * nx - K * nxy;
}
size_t Grid::j(size_t ijk)
{
    size_t K = k(ijk);
    return (ijk - K * nxy) / nx;
}
size_t Grid::k(size_t ijk)
{
    return ijk / nxy;
}
Coord Grid::ijk(const Coord &xyz)
{
    return Coord(i(xyz[1]), j(xyz[2]), k(xyz[3]));
}

// Domain Checks:
bool Grid::OutsideDomain(const Coord &xyz)
{
    return xyz[1] > endx or xyz[1] < startx or
           xyz[2] > endy or xyz[2] < starty or
           xyz[3] > endz or xyz[3] < startz;
}
bool Grid::OutsideDomain(double i, double j, double k)
{
    return i > nx - 1 or i < 0 or
           j > ny - 1 or j < 0 or
           k > nz - 1 or k < 0;
}

// Write Data to file:
void Grid::WriteFrametoCsv(float time, const RealBuffer &r, const RealBuffer &g, const RealBuffer &b, const RealBuffer &a, std::string directory, std::string name)
{
    PROFILE_FUNCTION();
    CreateDirectory(directory);

    name = (name == "") ? "data" : name;
    std::ofstream fileOut(directory + name + Format(time, 3) + "t " + std::to_string(nx) + "x " + std::to_string(ny) + "y" + ".csv");

    fileOut << "#nx=" << nx << "\n";
    fileOut << "#ny=" << ny << "\n";
    fileOut << "#nz=" << nz << "\n";
    fileOut << "#startx=" << startx << "\n";
    fileOut << "#starty=" << starty << "\n";
    fileOut << "#startz=" << startz << "\n";
    fileOut << "#endx=" << endx << "\n";
    fileOut << "#endy=" << endy << "\n";
    fileOut << "#endz=" << endz << "\n";
    fileOut << "#t=" << time << "\n";
    fileOut << "#x,y,z,r,g,b,a\n";

    for (size_t k = 0; k < nz; k++)
        for (size_t j = 0; j < ny; j++)
            for (size_t i = 0; i < nx; i++)
            {
                size_t ijk = Index(i, j, k);
                Coord x = xyz(i, j, k);
                fileOut << x[1] << "," << x[2] << "," << x[3] << ",";
                fileOut << r[ijk] << "," << g[ijk] << "," << b[ijk] << "," << a[ijk] << "\n";
            }

    fileOut.close();
}