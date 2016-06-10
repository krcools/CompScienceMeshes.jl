using Base.Test
using FixedSizeArrays
MUT = CompScienceMeshes

T = Float64
P = MUT.defaultpointtype(T,3)

width, height = 1.0, 1.0
mesh = MUT.meshrectangle(width, height, 0.5)

# This predicate tests whether a simplex is in the (x==0) plane
function pred(simplex)
    verts = mesh.vertices[simplex]
    for i in 1 : length(verts)
        if abs(verts[i][1]) > eps(Float64)
            return false
        end
    end
    return true
end


# Test if meshes can be intersected
line = MUT.meshsegment(width, 1/3, 3)
MUT.translate!(line, P(0.0, height, 0.0))
#line.vertices += Point(0.0, height, 0.0) # translation

edges = MUT.cells(mesh, 1)
γ1 = MUT.submesh(line, edges)
@test MUT.numcells(γ1) == 2
@test γ1.faces == [
    Vec(3,6),
    Vec(6,9)]

pred2 = MUT.interior_predicate(mesh)
γ2 = MUT.submesh(pred2, edges)
@test MUT.numcells(γ2) == 8

pred3 = MUT.overlap_predicate(mesh, line)
pred4 = x -> pred2(x) || pred3(x)
γ3 = MUT.submesh(pred4, edges)
@test MUT.numcells(γ3) == 10
