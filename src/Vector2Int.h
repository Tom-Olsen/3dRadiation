#ifndef __INCLUDE_GUARD_Edge_h__
#define __INCLUDE_GUARD_Edge_h__
#include <iostream>



struct Vector2Int
{
private:
    int data[2];
public:
    // Constructor:
    Vector2Int();
    Vector2Int(int x, int y);

    // Access:
    int& operator[](const int i);
    const int& operator[](const int i) const;

    // Basic Math:
    Vector2Int& operator=(const Vector2Int& e);

    // Boleans:
    bool operator==(const Vector2Int& other) const;
    bool operator!=(const Vector2Int& other) const;
    
    // Output:
    friend std::ostream& operator<<(std::ostream& os, const Vector2Int& e);
};
std::ostream& operator<<(std::ostream& os, const Vector2Int& e);
#endif //__INCLUDE_GUARD_Edge_h__