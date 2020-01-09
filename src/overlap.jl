import LinearAlgebra.cross
cross(a::Pt{2,T}, b::Pt{2,T}) where {T} = a[1]*b[2] - a[2]*b[1]

"""
Compute whether two flat patches of the same dimension overlap or not
"""
function overlap(p1::Simplex{U,D,C,N,T}, p2::Simplex{U,D,C,N,T}) where {U,D,C,N,T}

    # Are the patches in the same D-plane?
    throw(ErrorException("Not implemented yet!"))

end

###
# Version for points: overlap in this case is taken to mean that points approx conincide
###
function overlap(p::Simplex{U,0,C,1,T}, q::Simplex{U,0,C,1,T}) where {T,U,C}

    tol = sqrt(eps(T))
    return norm(p[1]-q[1]) < tol
end

"""
Compute whether two segments in 3D space overlap
"""
function overlap(p::Simplex{U,1,C,2,T}, q::Simplex{U,1,C,2,T}) where {T,U,C}

    tol = sqrt(eps(T))

    # Are these edges on the same line?
    o = q.vertices[2]
    u = p.vertices[1] - o
    t = q.tangents[1]
    t2 = dot(t,t)
    e = u - dot(u,t) / t2 * t

    # if vertex 1 of p is not in the line defined by q, they do not overlap
    norm(e) > tol && return false

    # if the tangents are not linearly dependent return false
    norm(cross(p.tangents[1],q.tangents[1])) > tol && return false

    # we determined p and q are on the same line. Do standard collision testing
    a1 = dot(p.vertices[1]-o, t)
    b1 = dot(p.vertices[2]-o, t)

    a1 < b1 ? (x1 = a1; y1 = b1) : (x1 = b1; y1 = a1)

    y1 <= zero(T) && return false # p to the left of q
    t2 <= x1 && return false      # q to the left of p

    return true
end

"""
Compute whether two triangles in 3D space overlap
"""
function overlap(p::Simplex{3,2,1,3,T}, q::Simplex{3,2,1,3,T}) where T

  tol = sqrt(eps(T))

  # Are the patches in the same plane?
  u1 = q.tangents[1]
  u2 = q.tangents[2]
  v = p.vertices[1] - q.vertices[2]

  # if vertex 1 of p is not in q, return false
  norm(dot(cross(u1,u2),v)) > tol && return false

  # if the two triangles are not coplanar return false
  norm(cross(p.normals[1], q.normals[1])) > tol && return false

  n = p.normals[1]
  for i in 1:3
    a = p.vertices[mod1(i+1,3)]
    b = p.vertices[mod1(i+2,3)]
    c = p.vertices[i]
    t = b - a
    m = cross(t,n)

    sp = zeros(T,3); sp[i] = dot(c-a, m)
    sq = T[dot(q.vertices[j]-a, m) for j in 1:3]

    minp, maxp = extrema(sp)
    minq, maxq = extrema(sq)

    maxq <= minp + tol && return false
    maxp <= minq + tol && return false
  end

  return true
end


function overlap(p::Simplex{3,3,0,4,T}, q::Simplex{3,3,0,4,T}) where T
    tol = sqrt(eps(T))
    for v in p.vertices
        u = carttobary(q,v)
        u = vcat(u, 1-sum(u))
        all(0+tol .<= u .<= 1-tol) && return true
    end
    u = carttobary(q, cartesian(center(p)))
    u = vcat(u, 1-sum(u))
    all(0+tol .<= u .<= 1-tol) && return true
    return false
end
