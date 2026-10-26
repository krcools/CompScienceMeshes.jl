using Test
using CompScienceMeshes
#Calling functions
refrect = meshrectangle(1.0, 1.0, 0.5, generator = :gmsh);
rect2 = meshrectangle(1.0, 1.0, 0.5, 2);
rect3 = meshrectangle(1.1, 0.6, 0.4, 3);

#Case: The function has a return
@test typeof(refrect) != Nothing
@test typeof(rect2) != Nothing
@test typeof(rect3) != Nothing

#Case: udim
@test length(rect2.vertices[1]) == 2
@test length(rect3.vertices[1]) == 3

#Case: udim = 3 is extension of udim = 2
@testset for i = 1:length(rect2.vertices)
    rect2.vertices[i][1] == rect3.vertices[i][1]
    rect2.vertices[i][2] == rect3.vertices[i][2]
end   
#@test rect2.faces == rect3.faces

#Case: The mesh returns vertices and faces
@test length(refrect.vertices) != 0
@test length(refrect.faces) != 0
@test length(rect2.vertices) != 0
@test length(rect2.faces) != 0
@test length(rect3.vertices) != 0
@test length(rect3.faces) != 0