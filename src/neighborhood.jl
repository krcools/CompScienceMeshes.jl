immutable MeshPointNM{U,D,C,N,T}
    patch::FlatCellNM{U,D,C,N,T}
    bary::SVector{D,T}
    cart::SVector{U,T}
end


Base.length(m::MeshPointNM) = length(m.cart)
Base.getindex(p::MeshPointNM, i::Int) = p.cart[i]

cartesian(mp::MeshPointNM) = mp.cart
parametric(mp::MeshPointNM) = mp.bary

"Return the barycentric coordinates of `mp`"
barycentric(mp::MeshPointNM) = SVector(mp.bary[1], mp.bary[2], 1-mp.bary[1]-mp.bary[2])

"""
A number defines a neighborhood in euclidian space
"""
jacobian(x::Number) = one(x)

jacobian(mp::MeshPointNM) = volume(mp.patch) * factorial(dimension(mp.patch))
tangents(mp::MeshPointNM, i) = mp.patch.tangents[i]
normal(mp::MeshPointNM) = mp.patch.normals[1]

function neighborhood(p::FlatCellNM, bary)
  D = dimension(p)
  T = coordtype(p)
  P = SVector{D,T}
  cart = barytocart(p, bary)
  MeshPointNM(p, P(bary), cart)
end

@generated function center{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T})
    uv = ones(T,D)/(D+1)
    :(neighborhood(p, $uv))
end
