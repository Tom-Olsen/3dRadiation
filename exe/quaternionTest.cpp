#include <iostream>
#include "../src/Stencil.h"
using namespace std;
using namespace glm;



vec3 NormalizedVec3(Tensor3 n)
{
    n = n.EuklNormalized();
    return vec3(n[1], n[2], n[3]);
}
void PrintVec3(const vec3& v, string name)
{ cout << name << "= (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << endl; };



void RotateStencil()
{
    vec3 from(0,0,1);
    vec3 to = NormalizedVec3(Tensor3(1,1,1));
    quat q = quat(from, to);
    quat qInv = quat(to, from);

    GaussLegendreStencil stencil(11);

    ofstream file0(OUTPUTDIR + "QuaternionRotationNormal.txt");
    ofstream file1(OUTPUTDIR + "QuaternionRotationRotated.txt");
    ofstream file2(OUTPUTDIR + "QuaternionRotationMixed.txt");
    file0 << "x, y, z, color" << endl;
    file1 << "x, y, z, color" << endl;
    file2 << "x, y, z, color" << endl;
    for(size_t d=0; d<stencil.nDir; d++)
    {
        Tensor3 p = stencil.Ct3(d);
        file0 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 0 << endl;
        p = q * stencil.Ct3(d);
        file1 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 1 << endl;
        p = glm::slerp(glm::quat(1,0,0,0),q,0.5f) * stencil.Ct3(d);
        file2 << p[1] << ", " << p[2] << ", " << p[3] << ", " << 2 << endl;
    }
    file0.close();
    file1.close();
    file2.close();

    cout << "3 files have been created in '" + OUTPUTDIR + "'. Plot them with ParaView (Filter:Table to Points)." << endl;
    cout << endl;
}



void QuaternionInversion()
{
    cout << "Quaternion Inversion:" << endl;
    vec3 from(0,0,1);
    vec3 to = NormalizedVec3(Tensor3(-1,1,1));
    quat q = quat(from, to);
    quat qInv = quat(to, from);

    vec3 v = NormalizedVec3(Tensor3(1,2,3));
    vec3 vRot = q * v;
    vec3 vInv0 = qInv * vRot;
    vec3 vInv1 = Invert(q) * vRot;

    PrintVec3(v    ,"    v");
    PrintVec3(vRot ," vRot");
    PrintVec3(vInv0,"vInv0");
    PrintVec3(vInv1,"vInv1");
    cout << endl;
}



void RecreateToFromQuatAndFrom()
{
    cout << "Recreate 'To' from 'Quat' and 'From':" << endl;
    vec3 from = NormalizedVec3(Tensor3(1,1,1));
    vec3 to = NormalizedVec3(Tensor3(0,1,1));
    quat q = quat(from,to);
    vec3 toRecreated = q * from;
    PrintVec3(from,"from");
    PrintVec3(to,"to");
    PrintVec3(toRecreated,"toRecreated");
    cout << endl;
}



void Slerp()
{
    cout << "Slerp:" << endl;
    vec3 from = NormalizedVec3(Tensor3(0,0,1));
    vec3 test = NormalizedVec3(Tensor3(1,1,0));
    vec3 to0 = NormalizedVec3(Tensor3(1,0,0));
    vec3 to1 = NormalizedVec3(Tensor3(0,1,0));
    quat q0 = quat(from,to0);
    quat q1 = quat(from,to1);

    int N = 15;
    quat q[N];

    for(int i=0; i<N; i++)
    {
        float t = i / (N - 1.0f);
        // q[i] = Normalized(slerp(q0,q1,t));
        q[i] = Normalized(q0 + (q1-q0) * t);
        qPrint(q[i], "q" + to_string(t));
    }
    cout << endl;

    for(int i=0; i<N; i++)
    {
        float t = i / (N - 1.0f);
        vec3 to = q[i] * from;
        PrintVec3(to,"to" + to_string(t));
    }
    cout << endl;

    cout << "x,y,z,color\n";
    cout << test[0] << "," << test[1] << "," << test[2] << "," << 1.5 << "\n";
    cout << 1 << "," << 0 << "," << 0 << "," << 1.5 << "\n";
    cout << 0 << "," << 1 << "," << 0 << "," << 1.5 << "\n";
    cout << 0 << "," << 0 << "," << 1 << "," << 1.5 << "\n";
    cout << 0 << "," << 0 << "," << 0 << "," << 1.5 << "\n";
    for(int i=1; i<N-1; i++)
    {
        float t = i / (N - 1.0f);
        vec3 to = q[i] * from;
        cout << to[0] << "," << to[1] << "," << to[2] << "," << t << "\n";
    }
    cout << endl;
}



int main()
{
    // RotateStencil();
    // QuaternionInversion();
    // RecreateToFromQuatAndFrom();
    Slerp();
}