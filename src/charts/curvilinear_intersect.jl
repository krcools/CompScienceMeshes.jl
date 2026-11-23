# Locating points in curvilinear triangles, and the overlap and intersection of two of them.

"""
    carttobary(ch::CurvilinearSimplex{U,2}, x; maxiter=20)

The parameters `u` of the point of `ch` closest to `x`: the solution of
`barytocart(ch, u) ≈ x` in the least squares sense by Gauss-Newton iteration, starting from the
parameters of `x` in the flat triangle through the corners of `ch`. For `x` on `ch` this is the
inverse of `barytocart`. The map of `ch` is a polynomial, defined beyond the reference triangle:
`u` lies outside of it when `x` lies on the continuation of `ch`.
"""
function carttobary(ch::CurvilinearSimplex{U,2,C,N,T}, x; maxiter=20) where {U,C,N,T}
    u = carttobary(simplex(nodes(ch)...), x)
    for _ in 1:maxiter
        t1, t2 = unitarybase(ch, u)
        r = barytocart(ch, u) - x
        G = SMatrix{2,2,T}(dot(t1, t1), dot(t2, t1), dot(t1, t2), dot(t2, t2))
        du = G \ SVector(dot(t1, r), dot(t2, r))
        u -= du
        norm(du) ≤ 4 * eps(T) && break
    end
    return u
end

# The parameters in `p` of the corners of `q` if `q` is an affine reparametrization of the map of
# `p` (continued beyond its reference triangle), as are a cell and its children in the barycentric
# refinement, or two children of the same cell; `nothing` otherwise. The corners determine the
# affine map, and the maps of `q` and of `p` composed with it are polynomials of the order of `q`
# that are equal if they agree at the interpolation nodes of `q`.
function _cornersinparams(p::CurvilinearSimplex{U,2,C,N,T}, q::CurvilinearSimplex{U,2,C,M,T,O}) where {U,C,N,T,M,O}
    us = map(c -> carttobary(p, c), nodes(q))
    scale = zero(T)
    for a in nodes(p), b in nodes(q)
        scale = max(scale, norm(a - b))
    end
    tol = sqrt(eps(T)) * scale
    for (r, t) in enumerate(gmsh_triangle_index_to_tuple[O])
        u = (t[1] * us[1] + t[2] * us[2] + t[3] * us[3]) / O
        norm(barytocart(p, u) - q.vertices[r]) ≤ tol || return nothing
    end
    return us
end

# The part of `p` over the triangle `(a, b, c)` of its parametric domain, as a chart of the same
# order, with the interpolation nodes of `p` mapped onto `(a, b, c)` affinely.
function _subchart(p::CurvilinearSimplex{U,2,C,N,T,O}, a, b, c) where {U,C,N,T,O}
    X = map(gmsh_triangle_index_to_tuple[O]) do t
        barytocart(p, (t[1] * a + t[2] * b + t[3] * c) / O)
    end
    return CurvilinearSimplex{U,2,C,N,T,O}(SVector{N}(X))
end

# The intersection of `q` with `p` as a polygon in the parametric domain of `p`, if `q` is an
# affine reparametrization of `p` (see `_cornersinparams`); empty otherwise.
function _intersectionparams(p::CurvilinearSimplex{U,2,C,N,T}, q::CurvilinearSimplex{U,2}) where {U,C,N,T}
    us = _cornersinparams(p, q)
    us === nothing && return SVector{2,T}[]
    reference = [SVector{2,T}(1, 0), SVector{2,T}(0, 1), SVector{2,T}(0, 0)]
    return sutherlandhodgman2d(collect(us), reference)
end

"""
    overlap(p::CurvilinearSimplex{U,2}, q::CurvilinearSimplex{U,2})

Whether the curvilinear triangles `p` and `q` share a part of positive area. This is decided
exactly when `q` is an affine reparametrization of the map of `p`, as for the cells of a mesh and
of its refinements; otherwise they are taken not to overlap, as non coplanar flat triangles.
"""
function overlap(p::CurvilinearSimplex{U,2,C,N,T}, q::CurvilinearSimplex{U,2}) where {U,C,N,T}
    poly = _intersectionparams(p, q)
    area = zero(T)
    for i in 3:length(poly)
        area += cross(poly[i-1] - poly[1], poly[i] - poly[1]) / 2
    end
    return abs(area) > 1e3 * eps(T)
end

"""
    intersection(p::CurvilinearSimplex{U,2}, q::CurvilinearSimplex{U,2})

The intersection of the curvilinear triangles `p` and `q` as curvilinear triangles of the order
of `p`, with its orientation, under the assumption of `overlap`: empty unless `q` is an affine
reparametrization of the map of `p`.
"""
function intersection(p::CurvilinearSimplex{U,2,C,N,T}, q::CurvilinearSimplex{U,2}) where {U,C,N,T}
    poly = _intersectionparams(p, q)
    pieces = typeof(p)[]
    for i in 3:length(poly)
        a, b, c = poly[1], poly[i-1], poly[i]
        area = cross(b - a, c - a) / 2
        abs(area) ≤ 1e3 * eps(T) && continue
        area < 0 && ((b, c) = (c, b))
        push!(pieces, _subchart(p, a, b, c))
    end
    return pieces
end
