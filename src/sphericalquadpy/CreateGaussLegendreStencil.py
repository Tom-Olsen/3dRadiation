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
    xyz = Q.xyz
    
    with open("../GaussLegendreStencil/GaussLegendreStencil" + str(order), "w") as f:
        for i in range(len(weights)):
            f.write("w[{}] =  {}; ".format(i, weights[i]))
            if xyz[i][0] >= 0:
                f.write("cx[{}] =  {}; ".format(i, xyz[i][0]))
            else:
                f.write("cx[{}] = {}; ".format(i, xyz[i][0]))
            if xyz[i][1] >= 0:
                f.write("cy[{}] =  {}; ".format(i, xyz[i][1]))
            else:
                f.write("cy[{}] = {}; ".format(i, xyz[i][1]))
            if xyz[i][2] >= 0:
                f.write("cz[{}] =  {};\n".format(i, xyz[i][2]))
            else:
                f.write("cz[{}] = {};\n".format(i, xyz[i][2]))

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

for order in range(3,20,2):
    GaußLegendreQuadratureToGrid(order)
    

# for i in range(1,21):
    # WriteGaussLegendreStencilToFile(i)