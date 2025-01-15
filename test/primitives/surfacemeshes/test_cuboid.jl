using Test
using CompScienceMeshes
#Calling functions
m_c = meshcuboid(1.0, 1.0, 1.0, 0.1);
mc= meshcuboid(1.0, 1.0, 1.0, 0.1, generator = :gmsh);
mcc = gmshcuboid(1.0, 1.0, 1.0, 0.1, physical = "OpenBox");

#Case: The function has a return
@test typeof(m_c) != Nothing
@test typeof(mc) != Nothing
@test typeof(mcc) != Nothing

#Case: The mesh returns vertices and faces
@test length(m_c.vertices) != 0
@test length(m_c.faces) != 0
@test length(mc.vertices) != 0
@test length(mc.faces) != 0
@test length(mcc.vertices) != 0
@test length(mcc.faces) != 0

#Case: Euler's formula for polyhedron
function calculate_edges(len, bdth, wid, h)
    n1 = Int(round(bdth/h))
    n2 = Int(round(len/h))
    n3 = Int(round(wid/h))
    n_e = (4 + 3*(n1 - 1) + 3*(n2 - 1) + 2*(n1 - 1)*(n2 - 1)) + 2*(
            3*n3 + 2*(n3*(n1 - 1))) + 2*(
            2*(n2 - 1)*n3 + n3) + (
            (n2 - 1)*(n1 - 1)*2 + (n2 - 1) + (n1 - 1)
            ) + 2*n2*n1 + 2*n3*n1 + 2*n2*n3
    return n_e
end
function euler(l, b, w, h)
    n = calculate_edges(l, b, w, h)
    v = length(meshcuboid(l, b, w, h).vertices)
    f = length(meshcuboid(l, b, w, h).faces)
    return v + f - n
end

@test euler(1.0, 1.0, 1.0, 0.01) == 2