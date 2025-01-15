using Test
using CompScienceMeshes

tt = tetmeshcuboid(1.0, 1.0, 1.0, 0.5);
t = tetmeshcuboid(1.0, 1.0, 1.0, 0.5, generator = :gmsh);

#Case: The function has a return
@test typeof(tt) != Nothing
@test typeof(t) != Nothing

#Case: The mesh returns vertices and faces
@test length(tt.vertices) != 0
@test length(tt.faces) != 0
@test length(t.vertices) != 0
@test length(t.faces) != 0