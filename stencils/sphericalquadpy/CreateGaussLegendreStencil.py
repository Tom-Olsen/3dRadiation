# source:
# https://github.com/thomascamminady/sphericalquadpy
import sphericalquadpy
import math
import statistics
import numpy as np



def CartesianToSpherical(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z / r)
    phi = math.atan2(y, x)
    if phi < 0:
        phi += 2 * math.pi
    return (r, theta, phi)

def ThetaListFromXyzList(xyz):
    theta = []
    for i in range(0, len(xyz)):
        x = xyz[i][0]
        y = xyz[i][1]
        z = xyz[i][2]
        r,th,ph = CartesianToSpherical(x,y,z)
        theta.append(th)
    return theta

def DistanceDeltaList(list):
    delta = []
    for i in range(0,len(list) - 1):
        delta.append(list[i+1] - list[i])
    return delta

def WriteGaussLegendreStencilToFile(order):
    Q = sphericalquadpy.gausslegendre.gausslegendre.GaussLegendre(order = order)
    weights = Q.weights
    weights = np.append(weights, 0)
    weights = np.append(weights, 0)
    
    xyz = Q.xyz
    xyz = np.append(xyz, [np.array([0, 0, 1])], axis=0)
    xyz = np.append(xyz, [np.array([0, 0, -1])], axis=0)
    
    theta = []
    phi = []
    for i in range(len(weights)):
        x = xyz[i][0]
        y = xyz[i][1]
        z = xyz[i][2]
        r,th,ph = CartesianToSpherical(x,y,z)
        theta.append(th)
        phi.append(ph)
        
    with open("../GaussLegendreStencil/GaussLegendreStencil" + str(order), "w") as f:
        f.write(f"this->nDir = {len(weights)} + this->nGhost;\n")
        f.write(f"AllocateBuffers();\n")
        for i in range(len(weights)):
            f.write(f"w[{i:>3}] = {weights[i]: .60f}; ")
            f.write(f"cx[{i:>3}] = {xyz[i][0]: .60f}; ")
            f.write(f"cy[{i:>3}] = {xyz[i][1]: .60f}; ")
            f.write(f"cz[{i:>3}] = {xyz[i][2]: .60f}; ")
            f.write(f"theta[{i:>3}] = {theta[i]: .60f}; ")
            f.write(f"phi[{i:>3}] = {phi[i]: .60f};\n")

def GaußLegendreQuadratureToGrid(order):
    Q = sphericalquadpy.gausslegendre.gausslegendre.GaussLegendre(order = order)

    # Get positions and weights:
    weights = Q.weights
    xyz = Q.xyz
    
    # Get unique theta list:
    thetaAll = ThetaListFromXyzList(xyz)
    theta = []
    for i in range(0,order):
        theta.append(thetaAll[2*i*order])
    theta.sort()

    # Distance lists:
    dTheta = DistanceDeltaList(theta)

    # Average Distance:
    dThetaAv= statistics.mean(dTheta)
    dThetaStDev = statistics.stdev(dTheta)

    # Print Stuff:
    # print("order:", order)
    # print("weights:", weights)
    # print("theta:", theta)
    # print("dTheta:", dTheta)
    print("dThetaAv:", dThetaAv)
    # print("dThetaStDev:", dThetaStDev)
    # print("Theta[0]", theta[0])
    # print("")

# for order in range(3,20,2):
    # GaußLegendreQuadratureToGrid(order)
    

for i in range(3,36,2):
    WriteGaussLegendreStencilToFile(i)
