#include "Grid.h"

// Constructors:
Grid::Grid(size_t nx_, size_t ny_, size_t nz_, Coord start_, Coord end_, size_t halo)
{
    if (end_[1] < start_[1])
        ExitOnError("Grid x€[a,b] has b<a!");
    if (end_[2] < start_[2])
        ExitOnError("Grid y€[a,b] has b<a!");
    if (end_[3] < start_[3])
        ExitOnError("Grid z[a,b] has b<a!");
    if (nx_ <= 1 or ny_ <= 1 or nz_ <= 1)
        ExitOnError("Grid must have at least 2 LP in each Dimension.");

    // Add 2 ghost cells and 1 extra cell due to off by one quirk.
    nx = nx_ + 2* halo;
    ny = ny_ + 2* halo;
    nz = nz_ + 2* halo;
    nxy = nx * ny;
    nxyz = nx * ny * nz;
    dx = (end_[1] - start_[1]) / (nx - 1.0 - 2.0 * halo);
    dy = (end_[2] - start_[2]) / (ny - 1.0 - 2.0 * halo);
    dz = (end_[3] - start_[3]) / (nz - 1.0 - 2.0 * halo);
    dt = m_cfl * std::min(std::min(dx, dy), dz);
    startx = start_[1] - dx;
    starty = start_[2] - dy;
    startz = start_[3] - dz;
    endx = end_[1] + dx;
    endy = end_[2] + dy;
    endz = end_[3] + dz;
}
Grid::Grid(const Grid &grid) :
nx(grid.nx), ny(grid.ny), nz(grid.nz), nxy(grid.nxy), nxyz(grid.nxyz),
dx(grid.dx), dy(grid.dy), dz(grid.dz),dt(grid.dt),
startx(grid.startx), starty(grid.starty), startz(grid.startz),
endx(grid.endx), endy(grid.endy), endz(grid.endz),
m_cfl(grid.m_cfl) {}

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
    return (x - startx) / dx;
}
double Grid::j(double y)
{
    return (y - starty) / dy;
}
double Grid::k(double z)
{
    return (z - startz) / dz;
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
    std::ofstream fileOut(directory + name + Format(time, 3) + "t " + std::to_string(nx) + "x " + std::to_string(ny) + "y " + std::to_string(nz) + "z.csv");

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

// Debugging:
void Grid::Print()
{
    std::cout << "Grid: nx=" << nx << "\n";
    std::cout << "Grid: ny=" << ny << "\n";
    std::cout << "Grid: nz=" << nz << "\n";
    std::cout << "Grid: nxy=" << nxy << "\n";
    std::cout << "Grid: nxyz=" << nxyz << "\n";
    std::cout << "Grid: dx=" << dx << "\n";
    std::cout << "Grid: dy=" << dy << "\n";
    std::cout << "Grid: dz=" << dz << "\n";
    std::cout << "Grid: dt=" << dt << "\n";
    std::cout << "Grid: startx=" << startx << "\n";
    std::cout << "Grid: starty=" << starty << "\n";
    std::cout << "Grid: startz=" << startz << "\n";
    std::cout << "Grid: endx=" << endx << "\n";
    std::cout << "Grid: endy=" << endy << "\n";
    std::cout << "Grid: endz=" << endz << "\n";
    std::cout << "Grid: cfl=" << m_cfl << "\n";
}