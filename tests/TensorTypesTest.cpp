#include <iostream>
#include "../src/TensorTypes.hh"

using namespace std;



int main()
{
    // Coord:
    Coord xyz(1.1,2.2,3.3);
    xyz.Print("xyz");
    PrintDouble(xyz[1],"x");
    PrintDouble(xyz[2],"y");
    PrintDouble(xyz[3],"z");
    cout << endl;

    // Tensor3:
    cout << "Tensor3" << endl;
    Tensor3 a3(1,2,3);
    a3.Print("a3");
    for(int i=1; i<4; i++)
        PrintDouble(a3[i],"a3[" + std::to_string(i) + "]");
    Tensor3 b3(4); 
    b3.Print("b3");
    for(int i=1; i<4; i++)
        PrintDouble(b3[i],"a3[" + std::to_string(i) + "]");
    cout << endl;
    
    // Tensor4:
    cout << "Tensor4" << endl;
    Tensor4 a4(1,2,3,4);
    a4.Print("a4");
    for(int i=0; i<4; i++)
        PrintDouble(a4[i],"a4[" + std::to_string(i) + "]");
    Tensor4 b4(4); 
    b4.Print("b4");
    for(int i=0; i<4; i++)
        PrintDouble(b4[i],"b4[" + std::to_string(i) + "]");
    cout << endl;

    // Tensor3x3:
    cout << "Tensor3x3" << endl;
    Tensor3x3 A33(0,1,2, 3,4,5, 6,7,8);
    A33.Print("A33");
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            PrintDouble(A33[{i,j}],"A33[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor3x3 B33(9);
    B33.Print("A33");
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            PrintDouble(B33[{i,j}],"A33[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;

    // Tensor4x4:
    cout << "Tensor4x4" << endl;
    Tensor4x4 A44(0,1,2,3, 2,3,4,5, 3,4,5,6, 4,5,6,7);
    A44.Print("A44");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            PrintDouble(A44[{i,j}],"A44[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor4x4 B44(1.234);
    B44.Print("A44");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            PrintDouble(B44[{i,j}],"A44[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;

    // Tensor3x3x3:
    cout << "Tensor3x3x3" << endl;
    Tensor3x3x3 A333(0);
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            for(int k=1; k<4; k++)
                A333[{i,j,k}] = (k-1) + 3*(j-1) + 9*(i-1);
    A333.Print("A333");
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            for(int k=1; k<4; k++)
                PrintDouble(A333[{i,j,k}],"A333[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;

    // Tensor4x4x4:
    cout << "Tensor4x4x4" << endl;
    Tensor4x4x4 A444(0);
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                A444[{i,j,k}] = k + 4*j + 16*i;
    A444.Print("A444");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                PrintDouble(A444[{i,j,k}],"A444[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;
}