export cartesian, jacobian, unormal, utangent, tangent, meshpoint, meshpoints
export MeshPoint, MeshPointNM

abstract MeshPoint{T}

cartesian(mp::MeshPoint) = mp.cart
jacobian(mp::MeshPoint) = volume(mp.cell) * factorial(dimension(mp.cell))
unormal(mp::MeshPoint) = mp.cell.unormal
utangent(mp::MeshPoint, i) = column(mp.cell.utangents, i)


immutable MeshPointNM{U,D,C,N,T} <: MeshPoint{T}
    patch::FlatCellNM{U,D,C,N,T}
    bary::Point{D,T}
    cart::Point{U,T}
end

tangent(mp::MeshPointNM, i) = mp.patch.tangents[i]
function utangent(mp::MeshPointNM, i)
    tang = mp.patch.tangents[i]
    return tang / norm(tang)
end
unormal(mp::MeshPointNM) = mp.patch.normals[1]
jacobian(mp::MeshPoint) = volume(mp.patch) * factorial(dimension(mp.patch))

function meshpoint{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}, bary)
    cart = barytocart(p, bary)
    #D = dimension(p)
    MeshPointNM(p, Point{D,T}(bary), cart)
end

function meshpoints{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}, uv::Array{T,2})
    numpoints = size(uv, 2)
    mps = Array(MeshPointNM{U,D,C,N,T}, numpoints)
    for i in 1:numpoints
        #bary = uv[:,i]
        mps[i] = meshpoint(p, uv[:,i])
    end
    return mps
end
