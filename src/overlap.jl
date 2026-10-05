import LinearAlgebra.cross
cross(a::Pt{2,T}, b::Pt{2,T}) where {T} = a[1]*b[2] - a[2]*b[1]

@inline _overlaptol(::Type{T}) where {T} = sqrt(eps(T))

@inline function _overlap_degenerate(length2, scale2)
    T = typeof(length2)
    return length2 <= eps(T) * max(length2, scale2)
end

"""
Compute whether two flat patches of the same dimension overlap or not
"""
function overlap(p1::Simplex{U,D,C,N,T}, p2::Simplex{U,D,C,N,T}) where {U,D,C,N,T}

    # Are the patches in the same D-plane?
    throw(ErrorException("Not implemented yet!"))

end

###
# For points, use the largest coordinate magnitude as the local scale. The unit
# floor keeps the tolerance useful for points close to the origin.
###
function overlap(p::Simplex{U,0,C,1,T}, q::Simplex{U,0,C,1,T}) where {T,U,C}
    scale = max(maximum(abs, p[1]), maximum(abs, q[1]), one(T))
    return norm(p[1] - q[1]) < _overlaptol(T) * scale
end

"""
Compute whether two segments in 3D space overlap

Segments that only share one endpoint are not considered overlapping; a
positive-length intersection is required.
"""
function overlap(p::Simplex{U,1,C,2,T}, q::Simplex{U,1,C,2,T}) where {T,U,C}

    tol = _overlaptol(T)

    # Are these edges on the same line?
    o = q.vertices[2]
    u = p.vertices[1] - o
    t = q.tangents[1]
    t2 = dot(t,t)
    p2 = dot(p.tangents[1], p.tangents[1])
    (_overlap_degenerate(t2, p2) || _overlap_degenerate(p2, t2)) && return false
    qlength = sqrt(t2)
    plength = sqrt(p2)
    length_scale = max(plength, qlength)
    e = u - dot(u,t) / t2 * t
    abstol = tol * length_scale

    # if vertex 1 of p is not in the line defined by q, they do not overlap
    norm(e) > abstol && return false

    # if the tangents are not linearly dependent return false
    norm(cross(p.tangents[1],q.tangents[1])) > tol * length_scale^2 && return false

    # we determined p and q are on the same line. Do standard collision testing
    a1 = dot(p.vertices[1]-o, t)
    b1 = dot(p.vertices[2]-o, t)

    a1 < b1 ? (x1 = a1; y1 = b1) : (x1 = b1; y1 = a1)

    intervaltol = tol * length_scale^2
    y1 - intervaltol <= zero(T) && return false # p to the left of q
    t2 <= x1 + intervaltol && return false      # q to the left of p

    return true
end

"""
Compute whether two triangles in 3D space overlap
"""
function overlap(p::Simplex{3,2,1,3,T}, q::Simplex{3,2,1,3,T}) where T

    tol = _overlaptol(T)
    edge_scale = zero(T)
    for i in 1:3, j in (i + 1):3
        edge_scale = max(edge_scale, norm(p.vertices[i] - p.vertices[j]))
        edge_scale = max(edge_scale, norm(q.vertices[i] - q.vertices[j]))
    end
    edge_scale2 = edge_scale * edge_scale
    _overlap_degenerate(edge_scale2, edge_scale2) && return false
    coplanartol = tol * edge_scale^3
    # `Simplex.normals` are unit normals, so their cross product is
    # dimensionless and must not be scaled by the triangle size.
    normaltol = tol
    # The separating-axis values below are dot products of an edge with a
    # cross product of an edge and a unit normal, and therefore scale like
    # length^2.
    separationtol = tol * edge_scale^2

  # Are the patches in the same plane?
  u1 = q.tangents[1]
  u2 = q.tangents[2]
  v = p.vertices[1] - q.vertices[2]

  # if vertex 1 of p is not in q, return false
    abs(dot(cross(u1,u2),v)) > coplanartol && return false

  # if the two triangles are not coplanar return false
    norm(cross(p.normals[1], q.normals[1])) > normaltol && return false

  n = p.normals[1]
  for i in 1:3
    a = p.vertices[mod1(i+1,3)]
    b = p.vertices[mod1(i+2,3)]
    c = p.vertices[i]
    t = b - a
    m = cross(t,n)

    pvalue = dot(c - a, m)
    qvalue1 = dot(q.vertices[1] - a, m)
    qvalue2 = dot(q.vertices[2] - a, m)
    qvalue3 = dot(q.vertices[3] - a, m)

    minp = min(zero(T), pvalue)
    maxp = max(zero(T), pvalue)
    minq = min(qvalue1, qvalue2, qvalue3)
    maxq = max(qvalue1, qvalue2, qvalue3)

    maxq <= minp + separationtol && return false
    maxp <= minq + separationtol && return false
  end

  return true
end


# Whether `x` lies in the closed tetrahedron `q`, shrunk by `tol` in barycentric
# coordinates. Barycentric coordinates are dimensionless, so `tol` needs no
# scaling with the size of the tetrahedron. Uses static arrays throughout.
@inline function _inside_tetrahedron(q::Simplex{3,3,0,4,T}, x, tol) where {T}
    t = q.tangents
    u = hcat(t[1], t[2], t[3]) \ (x - q.vertices[4])
    w = one(T) - (u[1] + u[2] + u[3])
    return tol <= u[1] <= 1 - tol && tol <= u[2] <= 1 - tol &&
           tol <= u[3] <= 1 - tol && tol <= w <= 1 - tol
end

"""
Compute whether two tetrahedra overlap: some vertex of `p`, or the center of `p`,
lies in the interior of `q`.
"""
function overlap(p::Simplex{3,3,0,4,T}, q::Simplex{3,3,0,4,T}) where T
    tol = _overlaptol(T)
    for v in p.vertices
        _inside_tetrahedron(q, v, tol) && return true
    end
    c = (p.vertices[1] + p.vertices[2] + p.vertices[3] + p.vertices[4]) / 4
    return _inside_tetrahedron(q, c, tol)
end



function overlap(p::Simplex{2,2,0,3,T}, q::Simplex{2,2,0,3,T}) where T

    #   tol = sqrt(eps(T))
    tol = 1e3 * eps(T)

  # Are the patches in the same plane?
  u1 = q.tangents[1]
  u2 = q.tangents[2]
  v = p.vertices[1] - q.vertices[2]

  for i in 1:3
    a = p.vertices[mod1(i+1,3)]
    b = p.vertices[mod1(i+2,3)]
    c = p.vertices[i]
    t = b - a
    m = StaticArrays.SVector{2,T}(t[2],-t[1])

    sp = zeros(T,3); sp[i] = dot(c-a, m)
    sq = T[dot(q.vertices[j]-a, m) for j in 1:3]

    minp, maxp = extrema(sp)
    minq, maxq = extrema(sq)

    maxq <= minp + tol && return false
    maxp <= minq + tol && return false
  end

  return true
end
