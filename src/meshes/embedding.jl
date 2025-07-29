"""
    embedding(small::AbstractMesh, big::AbstractMesh) -> S::AbstractSparseArray
    predicate(a,b) => extra constraints that must both be fullfiled by each
    cell `a` and `b` that are found to exist in `small` and `big` respectively.
    `a` and `b` are the indices of the corresponding edges. It is the user's
    responsibility to ensure that the predicate allows for all cells of `small`
    to be mapped to a cell in `big`. Specifying a predicate is useful for cases
    where cells in `small` matches/overlaps with multiple cells in `big`
    #TODO: change pred to accept simplices
"""
function embedding(small::AbstractMesh, big::AbstractMesh; predicate=(a,b)->true)

    @assert dimension(small) == 1
    @assert dimension(big)   == 1

    if length(big) == 0 || length(small) == 0
        return sparse(Int[], Int[], Float64[], length(small), length(big))
    end

    @assert length(big) > 0
    @assert length(small) > 0

    T = coordtype(small)
    tol = sqrt(eps(T))

    # store the edges of the big mesh in an octree
    ctrs_big = [cartesian(center(chart(big,v))) for v in big]
    tree_big = Octree(ctrs_big, zeros(size(ctrs_big)))
    
    cols = zeros(length(small))
    sgns = zeros(length(small))
    
    cells_big = collect(cells(big))
    for (i,c) in enumerate(small)
        V = chart(small, c).vertices
        ch = chart(small,c)
        x = cartesian(center(ch))
        pred(c,s) = fitsinbox(x, 0.0, c, s+tol)
        found = false
        for box in boxes(tree_big, pred)
            for j in box
                V′ = chart(big,j).vertices
                if (norm(V[1]-V′[1]) + norm(V[2]-V′[2]) < tol) && predicate(i,j)
                    cols[i] = j
                    sgns[i] = +1
                    found = true
                    break
                elseif (norm(V[1]-V′[2]) + norm(V[2]-V′[1]) < tol) && predicate(i,j)
                    cols[i] = j
                    sgns[i] = -1
                    found = true
                    break
                end
            end
            found && break
        end
    end

    # @assert !any(cols .== 0) "Unmapped cells exist"
    # @assert !any(sgns .== 0)

    I = findall(sgns .!= 0)
    rows = I
    cols = cols[I]
    sgns = sgns[I]
    # rows = range(1, length(small))[]

    # sparse(1:length(small), cols, sgns, length(small), length(big))
    return sparse(rows, cols, sgns, length(small), length(big))
end



function embedding(small::AbstractMesh{U,3}, big::AbstractMesh{U,3}; predicate=(a,b)->true) where{U}

    @assert dimension(small) == 2
    @assert dimension(big)   == 2

    if length(big) == 0 || length(small) == 0
        return sparse(Int[], Int[], Float64[], length(small), length(big))
    end

    @assert length(big) > 0
    @assert length(small) > 0

    T = coordtype(small)
    tol = sqrt(eps(T))

    # store the edges of the big mesh in an octree
    ctrs_big = [cartesian(center(chart(big,v))) for v in big]
    tree_big = Octree(ctrs_big, zeros(size(ctrs_big)))
    
    cols = zeros(length(small))
    sgns = zeros(length(small))
    
    cells_big = collect(cells(big))
    for (i,c) in enumerate(small)
        V = chart(small, c).vertices
        ch = chart(small,c)
        x = cartesian(center(ch))
        n = normal(center(ch))
        pred(c,s) = fitsinbox(x, 0.0, c, s+tol)
        found = false
        for box in boxes(tree_big, pred)
            for j in box
                ch′ = chart(big,j)
                V′ = ch′.vertices
                x′ = cartesian(center(ch′))
                n′ = normal(center(ch′))
                # if (norm(V[1]-V′[1]) + norm(V[2]-V′[2]) < tol) && predicate(i,j)
                norm(x-x′) >= tol && continue
                predicate(i,j) || continue
                if dot(n, n′) > 0
                    cols[i] = j
                    sgns[i] = +1
                    found = true
                    break
                # elseif (norm(V[1]-V′[2]) + norm(V[2]-V′[1]) < tol) && predicate(i,j)
                elseif dot(n, n′) < 0
                    cols[i] = j
                    sgns[i] = -1
                    found = true
                    break
                end
            end
            found && break
        end
    end

    # @assert !any(cols .== 0) "Unmapped cells exist"
    # @assert !any(sgns .== 0)

    I = findall(sgns .!= 0)
    rows = I
    cols = cols[I]
    sgns = sgns[I]
    # rows = range(1, length(small))[]

    # sparse(1:length(small), cols, sgns, length(small), length(big))
    return sparse(rows, cols, sgns, length(small), length(big))
end


@testitem "embedding for surfaces" begin

    m1 = meshrectangle(1.0, 1.0, 1.0, 3)
    m2 = CompScienceMeshes.translate(m1, point(1.0, 0.0, 0.0))

    m = CompScienceMeshes.weld(m1, m2)

    S = CompScienceMeshes.embedding(m1, m)

    # @show S
    @test Matrix(S) == [1 0 0 0; 0 1 0 0]
end


function embedding_topo(small::AbstractMesh, big::AbstractMesh)

    @assert dimension(small) == dimension(big) == 1
    @assert vertices(big) == vertices(small)

    @assert length(big) * length(small) > 0

    T = coordtype(small)
    tol = sqrt(eps(T))

    # store the edges of the big mesh in an octree
    ctrs_big = [cartesian(center(chart(big,v))) for v in big]
    tree_big = Octree(ctrs_big, zeros(size(ctrs_big)))
    
    # rows = 1:length(small)
    cols = zeros(length(small))
    sgns = zeros(length(small))
    
    cells_big = collect(cells(big))
    for (i,c) in enumerate(small)
        V = chart(small, c).vertices
        ch = chart(small,c)
        x = cartesian(center(ch))
        pred(c,s) = fitsinbox(x, 0.0, c, s+tol)
        found = false
        for box in boxes(tree_big, pred)
            for j in box
                c′ = cells_big[j]

                if c′ == c
                    cols[i] = j
                    sgns[i] = +1
                    found = true
                    break
                elseif c′ == reverse(c)
                    cols[i] = j
                    sgns[i] = -1
                    found = true
                    break
                end

                # V′ = chart(big,c′).vertices
                # ch′ = chart(big, c′)
                # if (norm(V[1]-V′[1]) + norm(V[2]-V′[2]) < tol) && predicate(i,j)
                #     cols[i] = j
                #     sgns[i] = +1
                #     found = true
                #     break
                # elseif (norm(V[1]-V′[2]) + norm(V[2]-V′[1]) < tol) && predicate(i,j)
                #     cols[i] = j
                #     sgns[i] = -1
                #     found = true
                #     break
                # end

            end
            found && break
        end
    end

    @assert !any(cols .== 0) "Unmapped cells exist"
    @assert !any(sgns .== 0)
    sparse(1:length(small), cols, sgns, length(small), length(big))
end
