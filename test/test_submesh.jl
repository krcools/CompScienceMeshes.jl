using Base.Test
using FixedSizeArrays
MUT = CompScienceMeshes

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
line.vertices += Point(0.0, height, 0.0)

edges = MUT.cells(mesh, 1)
γ = MUT.submesh(edges, line)
@test γ.faces == [
    Vec(3,6),
    Vec(6,9)]
