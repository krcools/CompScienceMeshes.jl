"""
permutate_mesh permutate the vertices of a mesh, while keeping the same cells.

Permutation is represented as a Vector v which sets the v[i]-th vertex at the i-th place.

mesh can be permutated from:
- Vector or SVector
- Permutations.Permutation
- another mesh with same universedimension
- another mesh with different universedimension, only works for combination of udims 2,3.

In the last two, control that the other mesh is a submesh of the starting mesh. In either direct


"""

function permutate(mesh::Mesh, σ::Union{Vector{<:Integer}, StaticArrays.SVector{Int32, <:Integer}})
    @assert numvertices(mesh) == length(σ)
    permutate(mesh, Permutations.Permutation(σ))
end

function permutate(mesh::Mesh, σ::Permutations.Permutation)
    mesh.vertices = mesh.vertices[σ.data]
    for (i, a) in enumerate(mesh.faces)
        mesh.faces[i] = σ'.data[a]
    end
end

function permutate(X::Mesh{U,D1}, Y::Mesh{U,D2}) where {U,D1,D2}
    tol = sqrt(eps(coordtype(X)))
    permut = Vector{Int32}()
    temp = collect(1:numvertices(X))
    for p in Y.vertices
        index = findfirst(isapprox(p;atol = tol), X.vertices)
        @assert !isnothing(index)

        temp[index] = 0
        append!(permut, index)
    end

    for i in temp
        if i !=0
            push!(permut, i)
    end end
    permutate(X, Permutations.Permutation(permut))
    X
end

function permutate_mesh(X::Mesh{2,D1}, Y::Mesh{3,D2}) where { D1, D2}
    tol = sqrt(eps(coordtype(X)))

    # assert that the z-coordinate for the vertices in X is constant.
    
    permut = Vector{Int32}()  
    temp = collect(1:numvertices(X))

    for p in Y.vertices
        q = zeros(eltype(X.vertices[begin]), 2)
        q[begin] = p[begin]
        q[2] = p[2]


        index = findfirst(isapprox(q;atol = tol), X.vertices)
        @assert !isnothing(index)   # vertex from Y exist in X
        @assert temp[index] != 0    # vertex from Y unique in X

        temp[index] = 0
        append!(permut, index)
    end

    for i in temp
        if i !=0
            push!(permut, i)
    end end
    
    permutate(X, Permutations.Permutation(permut))
    X
end

function permutate(X::Mesh{3,D1}, Y::Mesh{2,D2}) where {D1, D2}
    tol = sqrt(eps(coordtype(X)))

    # assert that the z-coordinate for the vertices in X is constant.
    for v in X.vertices
        @assert v[3] == X.vertices[1][3]
    end

    permut = Vector{Int32}()  
    temp = collect(1:numvertices(X))

    for p in Y.vertices
        q = zeros(eltype(X.vertices[begin]), 3)
        q[begin] = p[begin]
        q[2] = p[2]
        q[3] = X.vertices[1][3]

        index = findfirst(isapprox(q;atol = tol), X.vertices)
        @assert !isnothing(index)   # vertex from Y exist in X
        @assert temp[index] != 0    # vertex from Y unique in X

        temp[index] = 0
        append!(permut, index)
    end

    for i in temp
        if i !=0
            push!(permut, i)
    end end
    
    permutate_mesh(X, Permutations.Permutation(permut))
    X
end
