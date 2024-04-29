#include "Vector3Int.h"



// Constructor:
Vector3Int::Vector3Int()
{ data[0] = data[1] = data[2] = 0; }
Vector3Int::Vector3Int(int x, int y, int z)
{ data[0] = x; data[1] = y; data[2] = z; }

// Access:
int& Vector3Int::operator[](const int i)
{ return data[i]; }
const int& Vector3Int::operator[](const int i) const
{ return data[i]; }

// Basic Math:
Vector3Int& Vector3Int::operator=(const Vector3Int& t)
{
    data[0] = t[0];
    data[1] = t[1];
    data[2] = t[2];
    return *this;
}

// Booleans:
bool Vector3Int::operator==(const Vector3Int& other) const
{
    bool equalA = (data[0] == other[0] && data[1] == other[1] && data[2] == data[2]);
    bool equalB = (data[0] == other[1] && data[1] == other[2] && data[2] == data[0]);
    bool equalC = (data[0] == other[2] && data[1] == other[0] && data[2] == data[1]);
    return equalA || equalB || equalC;
}
bool Vector3Int::operator!=(const Vector3Int& other) const
{ return !(*this == other); }
    
// Output:
std::ostream& operator<<(std::ostream& os, const Vector3Int& t)
{
    os << t[0] << "," << t[1] << "," << t[2];
    return os;
}