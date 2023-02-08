# source:
# https://github.com/thomascamminady/sphericalquadpy
import csv
import sphericalquadpy
from numpy import exp 
import math

def f(x,y,z): return exp(-z**2)

def cartesian_to_spherical(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z / r)
    phi = math.atan2(y, x)
    return (r, theta, phi)

Q = sphericalquadpy.lebedev.Lebedev(order = 3)
# Q = sphericalquadpy.lebedev.Lebedev(order = 31)
# Q = sphericalquadpy.gausslegendre.gausslegendre.GaussLegendre(order = 13)
weights = Q.weights
xyz = Q.xyz

print(weights)
print(xyz)
# print(Q.integrate(f) - 9.38486877222642)

# xyz coordinates:
#with open('output.csv', mode='w') as csv_file:
#    writer = csv.writer(csv_file)
#    for i in range(len(weights)):
#        row = [weights[i], xyz[i][0], xyz[i][1], xyz[i][2]]
#        # row = [xyz[i][0], xyz[i][1], xyz[i][2]]
#        writer.writerow(row)
#
## theta,phi coorinates:
#with open('output.csv', mode='w') as csv_file:
#    writer = csv.writer(csv_file)
#    for i in range(len(weights)):
#        r, theta, phi = cartesian_to_spherical(xyz[i][0], xyz[i][1], xyz[i][2])
#        # row = [theta, phi]
#        row = [weights[i], theta, phi]
#        writer.writerow(row)