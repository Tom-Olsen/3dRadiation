#include <iostream>
#include "../src/Utility.hh"
using namespace std;



template<int order>
void MyAtan2()
{
    int n = 5;
    cout << "atan2f, MyAtan2, |delta|" << endl;
    for(int j=0; j<n; j++)
    {
        float y = -2.0 * j / (n-1.0) + 1.0;
        for(int i=0; i<n; i++)
        {
            float x = 2.0 * i / (n-1.0) - 1.0;
            cout << Format(atan2f(y,x)) << "," << Format(MyAtan2<order>(y,x)) << "," << Format(abs(atan2f(y,x) - MyAtan2<order>(y,x))) << endl;
        }
    }
    cout << endl;
}
template<int order>
void MySin()
{
    int n = 10;
    cout << "phi, sin(phi), MySin(phi), |delta|" << endl;
    for(int i=0; i<n; i++)
    {
        float phi = 2 * M_PI * i / (n - 1.0);
        cout << Format(phi) << "," << Format(sin(phi)) << "," << Format(MySin<order>(phi)) << "," << Format(sin(phi) - MySin<order>(phi)) << endl;
    }
    cout << endl;
}
template<int order>
void MyCos()
{
    int n = 10;
    cout << "phi, sin(phi), MySin(phi), |delta|" << endl;
    for(int i=0; i<n; i++)
    {
        float phi = 2 * M_PI * i / (n - 1.0);
        cout << Format(phi) << "," << Format(cos(phi)) << "," << Format(MyCos<order>(phi)) << "," << Format(cos(phi) - MyCos<order>(phi)) << endl;
    }
    cout << endl;
}



int main()
{
    MyAtan2<9>();
    MySin<9>();
    MyCos<9>();
}