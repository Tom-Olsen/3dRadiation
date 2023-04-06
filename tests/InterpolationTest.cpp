#include <iostream>
#include "../src/Grid.h"
#include "../src/Interpolation.hh"

using namespace std;



int main()
{
    size_t nx, ny, nz;
    nx = ny = nz = 50;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    RealBuffer data0(grid.nxyz);
    RealBuffer data1(grid.nxyz);
    RealBuffer data2(grid.nxyz);
    RealBuffer data3(grid.nxyz);
    for(size_t i=0; i<grid.nx; i++)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t k=0; k<grid.nz; k++)
    {
        size_t ijk = grid.Index(i,j,k);
        Coord xyz = grid.xyz(i,j,k);
        data0[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 1,0,0,0,0,0,0,1);
        data1[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,1,0,0,0,0,1,0);
        data2[ijk] = TrilinearInterpolation(i/49.0, j/49.0, k/49.0, 0,0,1,0,0,1,0,0);
        data3[ijk] = 1;
    }
    grid.WriteFrametoCsv(0,data0,data1,data2,data3,0,"output","Test_Interpolation");
}