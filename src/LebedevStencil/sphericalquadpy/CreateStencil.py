# source:
# https://github.com/thomascamminady/sphericalquadpy
import sphericalquadpy
from numpy import exp 


order = 15
Q = sphericalquadpy.lebedev.Lebedev(order = order)
weights = Q.weights
xyz = Q.xyz

print(weights)
print(xyz)

with open("../LebedevStencil" + str(order), "w") as f:
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