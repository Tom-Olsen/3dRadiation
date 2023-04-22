#include <iostream>
#include "../src/Grid.h"
using namespace std;



int main()
{
    size_t nx, ny, nz;
    nx = ny = nz = 20;
    Coord start(0,0,0);
    Coord end(1,1,1);
    Grid grid(nx, ny, nz, start, end);

    auto dist = [] (Coord a, Coord b)
    {
        double x = (a[1] - b[1]);
        double y = (a[2] - b[2]);
        double z = (a[3] - b[3]);
        return sqrt(x*x + y*y + z*z);
    };

    Coord rOrb(0.3,0.3,0.3);
    Coord gOrb(0.3,0.7,0.4);
    Coord bOrb(0.7,0.4,0.6);
    Coord aOrb(0.6,0.7,0.6);
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
        data0[ijk] = (dist(rOrb,xyz) < 0.2) ? 1 : 0;
        data1[ijk] = (dist(gOrb,xyz) < 0.2) ? 1 : 0;
        data2[ijk] = (dist(bOrb,xyz) < 0.2) ? 1 : 0;
        data3[ijk] = 1;
    }
    grid.WriteFrametoCsv (0,data0,data1,data2,data3,0,OUTPUTDIR,"GridTest");

    cout << "output/GridTest... has been created. Plot it with ParaView (Filter:Table to Points)." << endl;
}