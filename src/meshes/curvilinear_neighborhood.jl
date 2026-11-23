"""
    tangents(neighborhod, i) -> tangent_i

Return the i-th tangent vector at the neighborhood.
"""
function tangents(mp::MeshPointNM{T,C,D,U}, i) where {T,C<:CurvilinearSimplex,D,U}
    return unitarybase(mp.patch, mp.bary)[i]
end

# all tangents as the columns of a U×D matrix, as for flat simplices
function tangents(mp::MeshPointNM{T,C,D,U}) where {T,C<:CurvilinearSimplex,D,U}
    t = unitarybase(mp.patch, mp.bary)
    return hcat(t...)
end

function normal(mp::MeshPointNM{T,C,1,U}) where {T,C<:CurvilinearSimplex,U}
    t = tangents(mp, 1)
    nt = SVector(-t[2], t[1])

    return nt ./ norm(nt)
end

function normal(mp::MeshPointNM{T,C,2,U}) where {T,C<:CurvilinearSimplex,U}
    t1 = tangents(mp, 1)
    t2 = tangents(mp, 2)

    nt = cross(t1, t2)

    return nt ./ norm(nt)
end

function jacobian(mp::MeshPointNM{T,C,D,U}) where {T,C<:CurvilinearSimplex,D,U}
    return jacobian(mp.patch, mp.bary)
end
# The lazy neighborhood is meant for flat simplices, whose tangents and normal do not depend on
# the point and are taken from the chart. The tangents of a curved simplex do, and are computed
# from its unitary base at the point: as for quadrilaterals, the full neighborhood is returned.
neighborhood_lazy(ch::CurvilinearSimplex, u) = neighborhood(ch, u)
