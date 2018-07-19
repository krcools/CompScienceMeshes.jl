using CompScienceMeshes
using Test

m = meshrectangle(1.0,1.0,1.0,3)
p = chart(m, cells(m,i))

# bouding box of an array of points
c, s = boundingbox(p.vertices)
@test isa(c, eltype(p.vertices))
@test isa(s, Number)

# bouding box of a point
c, s = boundingbox(p.vertices[1])
@test c == p.vertices[1]
@test s == 0
