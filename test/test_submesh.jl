using Base.Test
MUT = CompScienceMeshes

mesh = MUT.meshrectangle(1.0, 1.0, 0.5)

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

sm = MUT.submesh(pred, mesh, 1)
#@test sm.faces == [1 2; 2 3]
@test sm.faces == [
    Vec(1,2),
    Vec(2,3)
]
