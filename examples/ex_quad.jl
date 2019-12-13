using CompScienceMeshes
using LinearAlgebra

o = point(3,4,2)
a = point(4,2,4)
b = point(2,8,2)
c = point(0,4,6)

tet = simplex(o,a,b,c)
qw = quadpoints(tet,4)
print(volume(tet),"==",abs(dot(a-o,cross(b-o,c-o)))/6,"==",abs(sum(w for (p,w) in qw)),"\n")

o,x,y,z = CompScienceMeshes.euclidianbasis(3)

tet = simplex(o,x,2y,3z)
volume(tet)

qw = quadpoints(tet,4)
length(qw)

print(abs(sum(w for (p,w) in qw)-volume(tet)) < 0.001,"\n")
print(sum(w*sin(p[1])*p[2]^2*p[3]^3 for (p,w) in qw))
#https://www.wolframalpha.com/input/?i=int_0^1+int_0^(2-2x)+int_0^(3-3x-1.5y)+sin(x)y^2z^3+dz+dy+dx
