
#export MeshPointNM

immutable MeshPointNM{U,D,C,N,T}
    patch::FlatCellNM{U,D,C,N,T}
    bary::Vec{D,T}
    cart::Vec{U,T}
end


paramtype{U,D,C,N,T}(::Type{MeshPointNM{U,D,C,N,T}}) = Vec{D,T}
pointtype{U,D,C,N,T}(::Type{MeshPointNM{U,D,C,N,T}}) = Vec{U,T}

Base.length(m::MeshPointNM) = length(m.cart)
Base.getindex(p::MeshPointNM, i::Int) = p.cart[i]

cartesian(mp::MeshPointNM) = mp.cart
parametric(mp::MeshPointNM) = mp.bary

"Return the barycentric coordinates of `mp`"
barycentric(mp::MeshPointNM) = Vec(mp.bary[1], mp.bary[2], 1-mp.bary[1]-mp.bary[2])

jacobian(mp::MeshPointNM) = volume(mp.patch) * factorial(dimension(mp.patch))
tangents(mp::MeshPointNM, i) = mp.patch.tangents[i]
function utangents(mp::MeshPointNM, i)
    tang = mp.patch.tangents[i]
    return tang / norm(tang)
end
normal(mp::MeshPointNM) = mp.patch.normals[1]

function meshpoint(p::FlatCellNM, bary)
  D = dimension(p)
  T = coordtype(p)
  P = Vec{D,T}
  cart = barytocart(p, bary)
  MeshPointNM(p, P(bary), cart)
end


meshpointtype{U,D,C,N,T}(::Type{FlatCellNM{U,D,C,N,T}}) =  MeshPointNM{U,D,C,N,T}
meshpointtype{U,D,C,N,T}(::FlatCellNM{U,D,C,N,T}) =  meshpointtype(FlatCellNM{U,D,C,N,T})


function meshpoints{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}, uv::Array{T,2})
    numpoints = size(uv, 2)
    mps = Array(MeshPointNM{U,D,C,N,T}, numpoints)
    for i in 1:numpoints
        mps[i] = meshpoint(p, uv[:,i])
    end
    return mps
end
