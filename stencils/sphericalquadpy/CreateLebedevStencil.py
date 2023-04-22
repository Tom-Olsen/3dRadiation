# source:
# https://github.com/thomascamminady/sphericalquadpy
import sphericalquadpy
import math


def CartesianToSpherical(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z / r)
    phi = math.atan2(y, x)
    if phi < 0:
        phi += 2 * math.pi
    return (r, theta, phi)

def WriteLebedevStencilToFile(order):
    Q = sphericalquadpy.lebedev.Lebedev(order = order)
    weights = Q.weights
    xyz = Q.xyz

    theta = []
    phi = []
    for i in range(len(weights)):
        x = xyz[i][0]
        y = xyz[i][1]
        z = xyz[i][2]
        r,th,ph = CartesianToSpherical(x,y,z)
        theta.append(th)
        phi.append(ph)

    with open("../LebedevStencil/LebedevStencil" + str(order), "w") as f:
        f.write(f"nDir = {len(weights)};\n")
        f.write(f"AllocateBuffers();\n")
        for i in range(len(weights)):
            f.write(f"w[{i:>3}] = {weights[i]: .60f}; ")
            f.write(f"cx[{i:>3}] = {xyz[i][0]: .60f}; ")
            f.write(f"cy[{i:>3}] = {xyz[i][1]: .60f}; ")
            f.write(f"cz[{i:>3}] = {xyz[i][2]: .60f}; ")
            f.write(f"theta[{i:>3}] = {theta[i]: .60f}; ")
            f.write(f"phi[{i:>3}] = {phi[i]: .60f};\n")
            
for i in range(3,32,2):
    WriteLebedevStencilToFile(i)
    
WriteLebedevStencilToFile(35)
WriteLebedevStencilToFile(41)
WriteLebedevStencilToFile(47)