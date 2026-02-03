using Test
using CompScienceMeshes
import Permutations

# check if permutated mesh has the same vertices and that the faces are the same.
function check_mesh(X, X_check)
    # test same vertices
    @test numvertices(X) == numvertices(X_check)
    @test Set(X.vertices) == Set(X_check.vertices)

    # test that each face has the same points in identical order
    @test numcells(X) == numcells(X_check)
    for i in 1:numcells(X)
        @test X.vertices[X.faces[i]] == X_check.vertices[X_check.faces[i]]
    end
end

# same universe dimension 
for u in 2:3
    local X = meshrectangle(1.0, 1.0, 0.5, u)
    local X_check = meshrectangle(1.0, 1.0, 0.5, u)
    local Y = meshrectangle(0.5, 0.5, 0.5, u)

    local σ = Permutations.Permutation(reverse(collect(1:numvertices(X))))

    permutate_mesh(X,σ.data)
    check_mesh(X, X_check)

    permutate_mesh(X, σ')
    @test X.vertices == X_check.vertices
    @test X.faces == X_check.faces

    permutate_mesh(X, Y)
    check_mesh(X, X_check)

    @test X.vertices[begin: begin+numvertices(Y)-1] == Y.vertices
end



# different universe dimension
let 
    local X = meshrectangle(1.0, 1.0, 0.5, 2)
    local X_check = meshrectangle(1.0, 1.0, 0.5, 2)
    local Y = meshrectangle(0.5, 0.5, 0.5, 3)

    permutate_mesh(X, Y)
    check_mesh(X, X_check)
    @test X.vertices[begin: begin+numvertices(Y)-1] == map(x->x[1:2], Y.vertices)


    permutate_mesh(X, X_check)

    translate!(Y, point(0,0,1))
    permutate_mesh(X, Y)
    check_mesh(X, X_check)
    @test X.vertices[begin: begin+numvertices(Y)-1] == map(x->x[1:2], Y.vertices)
end

let 
    local X = meshrectangle(1.0, 1.0, 0.5, 3)
    local X_check = meshrectangle(1.0, 1.0, 0.5, 3)
    local Y = meshrectangle(0.5, 0.5, 0.5, 2)
    translate!(X, point(0,0,1))
    translate!(X_check, point(0,0,1))

    permutate_mesh(X, Y)
    check_mesh(X, X_check)

    @test  map(x->x[1:2], X.vertices[begin: begin+numvertices(Y)-1]) == Y.vertices

end