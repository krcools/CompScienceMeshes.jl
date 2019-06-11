"""
    embedding(small::AbstractMesh, big::AbstractMesh) -> S::AbstractSparseArray
"""
function embedding(small::AbstractMesh, big::AbstractMesh)

    @assert dimension(small) == 1
    @assert dimension(big)   == 1

    T = coordtype(small)
    tol = sqrt(eps(T))

    # store the edges of the big mesh in an octree
    ctrs_big = [cartesian(center(chart(big,v))) for v in cells(big)]
    tree_big = Octree(ctrs_big, zeros(size(ctrs_big)))

    rows = 1:numcells(small)
    cols = zeros(numcells(small))
    sgns = zeros(numcells(small))

    cells_big = collect(cells(big))
    for (i,c) in enumerate(cells(small))
        # V = vertices(small,c)
        V = chart(small, c).vertices
        ch = chart(small,c)
        x = cartesian(center(ch))
        pred(c,s) = fitsinbox(x, 0.0, c, s+tol)
        found = false
        for box in boxes(tree_big, pred)
            for j in box
                c′ = cells_big[j]
                # V′ = vertices(big,c′)
                V′ = chart(big,c′).vertices
                ch′ = chart(big, c′)
                if norm(V[1]-V′[1]) + norm(V[2]-V′[2]) < tol
                    cols[i] = j
                    sgns[i] = +1
                    found = true
                    break
                elseif norm(V[1]-V′[2]) + norm(V[2]-V′[1]) < tol
                    cols[i] = j
                    sgns[i] = -1
                    found = true
                    break
                end
            end
            found && break
        end
    end

    @assert !any(cols .== 0)
    @assert !any(sgns .== 0)
    sparse(rows, cols, sgns, numcells(small), numcells(big))
end
