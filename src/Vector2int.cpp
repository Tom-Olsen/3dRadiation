#include "Vector2Int.h"



// Constructor:
Vector2Int::Vector2Int()
{ data[0] = data[1] = 0; }
Vector2Int::Vector2Int(int x, int y)
{ data[0] = x; data[1] = y; }

// Access:
int& Vector2Int::operator[](const int i)
{ return data[i]; }
const int& Vector2Int::operator[](const int i) const
{ return data[i]; }

// Basic Math:
Vector2Int& Vector2Int::operator=(const Vector2Int& e)
{
    data[0] = e[0];
    data[1] = e[1];
    return *this;
}

// Boleans:
bool Vector2Int::operator==(const Vector2Int& other) const
{
    bool equal = (data[0] == other[0] && data[1] == other[1]);
    bool inversEqual = (data[0] == other[1] && data[1] == other[0]);
    return equal || inversEqual;
}
bool Vector2Int::operator!=(const Vector2Int& other) const
{ return !(*this == other); }

// Output:
std::ostream& operator<<(std::ostream& os, const Vector2Int& e)
{
    os << e[0] << "," << e[1];
    return os;
}