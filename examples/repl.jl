using CompScienceMeshes
using FixedSizeArrays

p = Vec(1.0,0.0)
q = Vec(0.0,1.0)
r = Vec(-1.0000000000000002,-1.0)
a = Vec(1.0,0.0)
b = Vec(0.0,1.0)
c = Vec(0.0,0.0)

isct = sutherlandhodgman2d([p,q,r],[a,b,c])
