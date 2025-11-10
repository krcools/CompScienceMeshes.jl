"""
    tangents(neighborhod, i) -> tangent_i

Return the i-th tangent vector at the neighborhood.
"""
function tangents(mp::MeshPointNM{T,C,D,U}, i) where {T,C<:CurvilinearSimplex,D,U}
    return t = dpushforward(mp.patch, mp.bary[1])
end

function normal(mp::MeshPointNM{T,C,D,U}) where {T,C<:CurvilinearSimplex,D,U}
    t = tangents(mp, 1)
    nt = SVector(-t[2], t[1])

    return nt ./ norm(nt)
end

function jacobian(mp::MeshPointNM{T,C,D,U}) where {T,C<:CurvilinearSimplex,D,U}
    return jacobian(mp.patch, mp.bary[1])
end