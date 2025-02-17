using Test
using CompScienceMeshes
radius = 1.0;
# md2 = meshsphere(radius, 0.1, delaunay =:(2D));
# md3 = meshsphere(radius, 0.1, delaunay =:(3D));
mg = meshsphere(radius, 0.1, generator = :gmsh);

#Case: The function has a return
# @test typeof(md2) != Nothing
# @test typeof(md3) != Nothing
@test typeof(mg) != Nothing

#Case: The mesh returns vertices and faces
# @test length(md2.vertices) != 0
# @test length(md2.faces) != 0
# @test length(md3.vertices) != 0
# @test length(md3.faces) != 0
@test length(mg.vertices) != 0
@test length(mg.faces) != 0
