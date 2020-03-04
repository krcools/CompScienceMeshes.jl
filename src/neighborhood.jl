struct MeshPointNM{U,D,C,N,T}
    patch::Simplex{U,D,C,N,T}
    bary::SVector{D,T}
    cart::SVector{U,T}
end


Base.length(m::MeshPointNM) = length(m.cart)
Base.getindex(p::MeshPointNM, i::Int) = p.cart[i]

"""
    cartesian(neighborhood) -> point

Return the cartesian coordinates of the point where `neighborhood` is located.
"""
function cartesian(mp::MeshPointNM)
    mp.cart
end

"""
    parametric(neighborhood) -> point

Return the parameters where `neighborhood` is located.
"""
function parametric(mp::MeshPointNM)
    mp.bary
end

chart(nbd::MeshPointNM) = mp.patch

"Return the barycentric coordinates of `mp`"
barycentric(mp::MeshPointNM) = SVector(mp.bary[1], mp.bary[2], 1-mp.bary[1]-mp.bary[2])

"""
A number defines a neighborhood in euclidian space
"""
jacobian(x::Number) = one(x)

jacobian(mp::MeshPointNM) = volume(mp.patch) * factorial(dimension(mp.patch))

"""
    tangents(neighborhod, i) -> tangent_i

Return the i-th tangent vector at the neighborhood.
"""
function tangents(mp::MeshPointNM, i)
    mp.patch.tangents[i]
end


"""
    normal(neighborhood)

Return the normal at a neighborhood on a surface.
"""
function normal(mp::MeshPointNM)
    mp.patch.normals[1]
end

"""
    neighborhood(chart, params)

Create a neighborhood from a chart and a set of parameter values.
"""
function neighborhood(p::Simplex, bary)
  D = dimension(p)
  T = coordtype(p)
  P = SVector{D,T}
  cart = barytocart(p, bary)
  MeshPointNM(p, P(bary), cart)
end

"""
    center(simplex) -> neighborhood

Create the neighborhood at the center of the simplex, i.e. the point corresponding
to parameter `(1/(D+1), 1/(D+1), ...)` where `D` is the simplex dimension.
"""
@generated function center(p::Simplex{U,D,C,N,T}) where {U,D,C,N,T}
    uv = ones(T,D)/(D+1)
    :(neighborhood(p, $uv))
end
