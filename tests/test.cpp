#include <iostream>
#include <span>
#include <vector>
#include "../src/TensorTypes.hh"

using namespace std;



void A()
{
    std::vector<int,AlignedArrayAllocator<int>> vec = {1,2,3,4,5,6,7,8,9,10};

    for(auto it = begin(vec); it != end(vec); it++)
        cout << *it << endl;

    span<int> subVec = std::span(&vec[4], std::next(&vec[7]));
    
    for(auto it = begin(subVec); it != end(subVec); it++)
        cout << *it << endl;
}
void B()
{
    UnstructuredMatrix<int> mat;
    vector<int> a = {1,2,3};
    vector<int> b = {10,20};
    vector<int> c = {5,6,7,8,9};

    mat.AddRow(a);
    mat.AddRow(b);
    mat.AddRow(c);

    mat.Print();
}

int main()
{
    A();
    B();
}