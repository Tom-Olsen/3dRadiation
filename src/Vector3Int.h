#ifndef __INCLUDE_GUARD_Vector3Int_h__
#define __INCLUDE_GUARD_Vector3Int_h__
#include <iostream>
#include "Vector3.h"



struct Vector3Int
{
private:
    int data[3];
public:
    // Constructor:
    Vector3Int();
    Vector3Int(int x, int y, int z);

    // Access:
    int& operator[](const int i);
    const int& operator[](const int i) const;
    
    // Basic Math:
    Vector3Int& operator=(const Vector3Int& t);
    
    // Boleans:
    bool operator==(const Vector3Int& other) const;
    bool operator!=(const Vector3Int& other) const;

    // Output:
    friend std::ostream& operator<<(std::ostream& os, const Vector3Int& t);
};
std::ostream& operator<<(std::ostream& os, const Vector3Int& t);
#endif //__INCLUDE_GUARD_Vector3Int_h__