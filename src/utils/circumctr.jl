# import CompScienceMeshes.Simplex

# circumcenter(splx) = circumcenter(splx, dimtype(splx))

# function circumcenter(simplex, dt)

#     F = faces(simplex)
#     c1 = circumcenter(F[1], dimtype(F[1]))
#     c2 = circumcenter(F[2], dimtype(F[2]))

#     n1 = normal(F[1])
#     n2 = normal(F[2])
#     N = [n1 -n2]

#     T = tangents(simplex)

#     c = -c1 + c2
#     s = (T'*N) \ (T'*c)
#     return c1 + n1 * s[1] 
# end

# function circumcenter(simplex, ::Type{Val{1}})
#     return (simplex.vertices[1] + simplex.vertices[2])/2
# end

function circumcenter(simplex)

    V = vertices(simplex)
    G = V'*V
    T = coordtype(simplex)
    L = @SMatrix ones(T, size(G,1), 1)
    Z = @SMatrix zeros(T,1,1)

    A = [2G L; L' Z]


    D = vec(sum(V.^2, dims=1))
    z = @SVector zeros(T,1)

    b = [D; z]
    x = A \ b

    c = V*x[1:end-1]

    return c
end