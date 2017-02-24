using Base.Test

using CompScienceMeshes
using StaticArrays

T = Float64
#P = CompScienceMeshes.defaultpointtype(T,3)
P = SVector{3,T}

width, height = 1.0, 1.0
mesh = meshrectangle(width, height, 0.5)

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
line = meshsegment(width, 1/3, 3)
translate!(line, P(0.0, height, 0.0))

edges = skeleton(mesh, 1)
γ1 = submesh(line, edges)
@test numcells(γ1) == 2
@test γ1.faces == [
    index(3,6),
    index(6,9)]

pred2 = interior_predicate(mesh)
γ2 = submesh(pred2, edges)
@test numcells(γ2) == 8

pred3 = overlap_predicate(mesh, line)
pred4 = x -> pred2(x) || pred3(x)
γ3 = submesh(pred4, edges)
@test numcells(γ3) == 10
