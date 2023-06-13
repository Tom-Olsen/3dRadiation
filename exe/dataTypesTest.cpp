#include <iostream>
#include "../src/DataTypes.hh"
using namespace std;



void TestCoord()
{
    Coord xyz(1,1,1);
    xyz.Print("xyz");
    cout << "xyz:" << xyz << endl;
    PrintDouble(xyz[1],"x");
    PrintDouble(xyz[2],"y");
    PrintDouble(xyz[3],"z");
    cout << endl;
}
void TestTensor3()
{
    cout << "Tensor3" << endl;
    Tensor3 a(1,2,3);
    a.Print("a3");
    cout << "a:" << a << endl;;
    for(int i=1; i<4; i++)
        PrintDouble(a[i],"a[" + std::to_string(i) + "]");
    Tensor3 b(4); 
    b.Print("b");
    cout << "b:" << b << endl;;
    for(int i=1; i<4; i++)
        PrintDouble(b[i],"b[" + std::to_string(i) + "]");
    cout << endl;
}
void TestTensor4()
{
    cout << "Tensor4" << endl;
    Tensor4 a(1,2,3,4);
    a.Print("a");
    cout << "a:" << a << endl;;
    for(int i=0; i<4; i++)
        PrintDouble(a[i],"a[" + std::to_string(i) + "]");
    Tensor4 b(4); 
    b.Print("b");
    cout << "b:" << b << endl;;
    for(int i=0; i<4; i++)
        PrintDouble(b[i],"b[" + std::to_string(i) + "]");
    cout << endl;
}
void TestTensor3x3()
{
    cout << "Tensor3x3" << endl;
    Tensor3x3 A(0,1,2, 3,4,5, 6,7,8);
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            PrintDouble(A[{i,j}],"A[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor3x3 B(9);
    B.Print("B");
    cout << "B:" << B << endl;;
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            PrintDouble(B[{i,j}],"B[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;
}
void TestTensor4x4()
{
    cout << "Tensor4x4" << endl;
    Tensor4x4 A(0,1,2,3, 2,3,4,5, 3,4,5,6, 4,5,6,7);
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            PrintDouble(A[{i,j}],"A[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor4x4 B(1.234);
    B.Print("B");
    cout << "B:" << B << endl;;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            PrintDouble(B[{i,j}],"B[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;
}
void TestTensor3x3x3()
{
    cout << "Tensor3x3x3" << endl;
    Tensor3x3x3 A(0);
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            for(int k=1; k<4; k++)
                A[{i,j,k}] = (k-1) + 3*(j-1) + 9*(i-1);
    A.Print("A");
    cout << "A:" << A << endl;
    for(int i=1; i<4; i++)
        for(int j=1; j<4; j++)
            for(int k=1; k<4; k++)
                PrintDouble(A[{i,j,k}],"A[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;
}
void TestTensor4x4x4()
{
    cout << "Tensor4x4x4" << endl;
    Tensor4x4x4 A(0);
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                A[{i,j,k}] = k + 4*j + 16*i;
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                PrintDouble(A[{i,j,k}],"A[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;
}



void TestUnstructuredMatrix()
{
    UnstructuredMatrix<double> unstMat;
    double a[3] = {0,1,2};
    double b[4] = {10,11,12,13};
    double c[2] = {20,21};
    double d[5] = {30,31,32,33,34};

    unstMat.AddRow(a,3);
    unstMat.AddRow(b,4);
    unstMat.AddRow(c,2);
    unstMat.AddRow(d,5);

    unstMat.Print();
    unstMat.PrintFlat();

    for(int d=0; d<4; d++)
    {
        std::span row = unstMat.Row(d);
        cout << d << ": ";
        for(auto it = begin(row); it != end(row); it++)
            cout << *it << ",";
        cout << endl;
    }

    cout << "rowCount = " << unstMat.RowCount() << endl;
}



int main()
{
    TestCoord();
    TestTensor3();
    TestTensor4();
    TestTensor3x3();
    TestTensor4x4();
    TestTensor3x3x3();
    TestTensor4x4x4();
    TestUnstructuredMatrix();
}