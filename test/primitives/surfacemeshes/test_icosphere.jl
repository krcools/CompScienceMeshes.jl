using Test
using CompScienceMeshes
#Calling functions
ico1 = meshicosphere(1);
ico2 = meshicosphere(2);
ico3 = meshicosphere(2, 2.0);
ico4 = meshicosphere(2, 1.0, 30);

#Case: Less than or equal to the number of required vertices
@test ico2.vertices == ico4.vertices
@test ico2.faces == ico4.faces

#Case: Scaling 
@test length(ico2.vertices) == length(ico3.vertices)
@test length(ico2.faces) == length(ico3.faces)
@test ico2.vertices != ico3.vertices
@test ico2.faces == ico3.faces

#Case: The function has a return
@test typeof(ico1) != Nothing
@test typeof(ico2) != Nothing
@test typeof(ico3) != Nothing
@test typeof(ico4) != Nothing

#Case: Comparison on the sizes of mesh
@test length(ico1.vertices) < length(ico2.vertices)

#Case: The mesh returns vertices and faces
@test length(ico1.vertices) != 0
@test length(ico1.faces) != 0