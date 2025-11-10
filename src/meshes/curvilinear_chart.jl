#TODO: Extend to 2D and 3D simplices.

"""
    CurvilinearSimplex{U,D,C,N,T}

A curvilinear simplex (currently specialized to 1D) embedded in `U`-dimensional
space, represented by `N` interpolation vertices and `N` parametric nodes.

# Fields
- `vertices :: SVector{N, SVector{U,T}}`  
    Coordinates of the interpolation/control points in physical space.

- `ζnodes  :: SVector{N,T}`
    Parametric locations of the nodes on the reference interval `[0,1]` used
    for high-order interpolation (e.g. Gmsh’s high-order node ordering).

# Description
A `CurvilinearSimplex` represents a single high-order 1D finite-/boundary-element.
Its geometry is given by an isoparametric map

x(ζ) = Σ_r ℓ_r(ζ) * vertices[r],

where `ℓ_r(·)` are the Lagrange polynomials defined by `ζnodes`.

The type parameters are:
- `U`: ambient (universe) dimension,
- `D`: topological dimension (here always `1`),
- `C`: codimension (= `U-D`),
- `N`: number of local nodes,
- `T`: coordinate type (e.g. `Float64`).

This is the geometric building block used for curvilinear 1D meshes.
"""
struct CurvilinearSimplex{U,D,C,N,T} <: AbstractSimplex{U,D}
    vertices::SVector{N,SVector{U,T}}
    ζnodes::SVector{N,T}
end

"""
    nodes(s::CurvilinearSimplex{U,1})

Return the physical coordinates of the two endpoints of a 1D curvilinear element,
identified by the parametric values `ζ = 1` (left) and `ζ = 0` (right) under the
barycentric convention `ζ = λ_left`.

# Returns
`(x_left, x_right) :: (SVector{U,T}, SVector{U,T})`
"""
@inline function nodes(s::CurvilinearSimplex{U,1}) where {U}
    return SVector(s.vertices[begin], s.vertices[begin + 1])
end

# convenience
vertices(s::CurvilinearSimplex) = s.vertices

"""
    coordtype(s::CurvilinearSimplex)
    coordtype(::Type{CurvilinearSimplex})

Return the scalar coordinate type `T` used to represent vertex coordinates.
"""
coordtype(::CurvilinearSimplex{U,D,C,N,T}) where {U,D,C,N,T} = T

"""
    dimension(s::CurvilinearSimplex)
    dimension(::Type{CurvilinearSimplex})

Return the topological dimension of the simplex (`1` for line elements).
"""
dimension(::CurvilinearSimplex{U,D,C,N,T})  where {U,D,C,N,T} = D

"""
    universedimension(s::CurvilinearSimplex)
    universedimension(::Type{CurvilinearSimplex})

Return the dimension of the embedding space (e.g. `2` for curves in the plane).
"""
universedimension(::CurvilinearSimplex{U,D,C,N,T}) where {U,D,C,N,T} = U

"""
    vertextype(s::CurvilinearSimplex)
    vertextype(::Type{CurvilinearSimplex})

Return the static vector type used to represent vertex coordinates.
"""
vertextype(::CurvilinearSimplex{U,D,C,N,T}) where {U,D,C,N,T} = SVector{U,T}

# type-form
coordtype(::Type{CurvilinearSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = T
dimension(::Type{CurvilinearSimplex{U,D,C,N,T}})  where {U,D,C,N,T} = D
universedimension(::Type{CurvilinearSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = U
vertextype(::Type{CurvilinearSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = SVector{U,T}

"""
    neighborhood(ch::CurvilinearSimplex, bary)

Construct a `MeshPoint` located at the parametric (barycentric) coordinate `bary`
on the given curvilinear simplex.

The returned object stores:
- the chart `ch`,
- the parametric point `bary`,
- the corresponding physical coordinates `cartesian`.

This is the canonical way to create a point for evaluating reference functions
or geometric quantities on curvilinear elements.
"""
function neighborhood(p::C, bary) where {C<:CurvilinearSimplex}
    D = dimension(p)
    T = coordtype(p)
    P = SVector{D,T}
    cart = barytocart(p, T.(bary))
    Q = typeof(cart)
    U = length(cart)
    MeshPointNM{T,C,D,U}(p, P(bary), cart)
end

barytocart(ch::CurvilinearSimplex{U,1}, u) where {U} = pushforward(ch, u[1])

cartesian(ch::AbstractSimplex{U,1}, u::SVector{1,<:Real})  where {U} = barytocart(ch, u)
parametric(ch::AbstractSimplex{U,1}, u::SVector{1,<:Real}) where {U} = u

"""
    pushforward(ch::CurvilinearSimplex{U,1,C,N,T}, ζ)

Evaluate the isoparametric mapping from the scalar barycentric coordinate
`ζ = λ_left ∈ [0,1]` to physical space:

x(ζ) = Σ_r ℓ_r(ζ) * ch.vertices[r],

where `ℓ_r` are the Lagrange polynomials associated with `ch.ζnodes`.

With this convention, `ζ=1` maps to the left endpoint and `ζ=0` to the right.
"""
@inline function pushforward(ch::CurvilinearSimplex{U,1,C,N,T}, ζ::Real) where {U,C,N,T}
    z = T(ζ)
    s = zero(SVector{U,T})
    @inbounds for r in 1:N
        s += ch.vertices[r] * _ℓ(ch.ζnodes, r, z)
    end
    s
end

"""
    dpushforward(ch::CurvilinearSimplex, ζ)

Return the derivative of the isoparametric map `x(ζ)` with respect to `ζ`.

Mathematically:

dx/dζ = Σ_r ℓ'_r(ζ) * ch.vertices[r].


This quantity is used for computing Jacobians, tangent vectors, and element
volumes on curvilinear 1D elements.
"""
@inline function dpushforward(ch::CurvilinearSimplex{U,1,C,N,T}, ζ::Real) where {U,C,N,T}
    z = T(ζ)
    s = zero(SVector{U,T})
    @inbounds for r in 1:N
        s += ch.vertices[r] * _dℓ(ch.ζnodes, r, z)
    end
    s
end

"""
    _ℓ(ζnodes, r, ζ)

Evaluate the `r`-th 1D Lagrange basis polynomial `ℓ_r(ζ)` associated with the
nodes given in `ζnodes`.

The polynomial satisfies:
- `ℓ_r(ζnodes[r]) = 1`
- `ℓ_r(ζnodes[j]) = 0` for all `j ≠ r`

This is the standard explicit Lagrange form:

ℓ_r(ζ) = ∏_{j ≠ r} (ζ - ζ_j) / (ζ_r - ζ_j).

# Internal
Used for geometric mapping, Jacobians, and interpolation on curvilinear 1D elements.
"""
@inline function _ℓ(ζnodes::SVector{D,T}, r::Int, ζ::T) where {D,T<:Real}
    ζr = ζnodes[r]
    num = one(T); den = one(T)
    @inbounds for j in 1:D
        j == r && continue
        num *= (ζ - ζnodes[j])
        den *= (ζr - ζnodes[j])
    end

    return num/den
end

"""
    _dℓ(ζnodes, r, ζ)

Evaluate the derivative of the `r`-th Lagrange basis polynomial at ζ.
This is the analytical derivative of `_ℓ`.

Mathematically:

ℓ'r(ζ) = Σ{k ≠ r}
∏{j ≠ r,k} (ζ - ζ_j)
/ (ζ_r - ζ_k) / ∏{j ≠ r} (ζ_r - ζ_j).


This is used to compute the derivative of the physical mapping and hence the
Jacobian and tangents.
"""
@inline function _dℓ(ζnodes::SVector{D,T}, r::Int, ζ::T) where {D,T<:Real}
    ζr = ζnodes[r]
    s = zero(T)
    @inbounds for k in 1:D
        k == r && continue
        num = one(T); den = (ζr - ζnodes[k])
        @inbounds for j in 1:D
            (j == r || j == k) && continue
            num *= (ζ - ζnodes[j])
            den *= (ζr - ζnodes[j])
        end
        s += num/den
    end
    s
end


"""
    gmsh_line_ζnodes(::Val{N})

Return the parametric node locations `ζ ∈ [0,1]` for a 1D high-order line element,
arranged to match the element’s local node numbering (compatible with Gmsh element
node indices and the barycentric convention used here).

By convention in this code, the scalar parameter is the *left* barycentric
coordinate: `ζ = 1` corresponds to the left endpoint and `ζ = 0` to the right
endpoint.

- For `N = 2` (linear): returns `(1.0, 0.0)` → `(left, right)`.
- For `N > 2`: returns `(1.0, 0.0, ζ₃, …, ζ_N)` where the interior nodes `ζ_k`
  lie strictly between 0 and 1. (Their order is chosen to be consistent with
  Gmsh’s line-element node numbering and this code’s interpolation routines.)

This tuple is used by the Lagrange map `x(ζ) = Σ_r ℓ_r(ζ) * vertices[r]`.
"""
gmsh_line_ζnodes(::Val{N}) where {N} =
    N == 2 ? SVector{2,Float64}(1.0, 0.0) :
             SVector{N,Float64}(1.0, 0.0, (k/(N-1) for k in reverse(1:N-2))...)


function chart(
    m::CompScienceMeshes.CurvilinearMesh{U,N,T,O},
    faceid::Int
) where {U,N,T,O}

    vindices = m.faces[faceid]
    X = m.vertices[vindices]
    ζnodes = gmsh_line_ζnodes(Val(N))                   # match that order

    CompScienceMeshes.CurvilinearSimplex{U,1,U-1,N,T}(X, ζnodes)
end

"""
    jacobian(ch, ζ)

Return the scalar Jacobian `|dx/dζ|` of the 1D isoparametric mapping at
parametric coordinate `ζ`.
"""
jacobian(ch::CurvilinearSimplex{U,1}, ζ::Real) where {U} = norm(dpushforward(ch, ζ))

"""
    center(ch::CurvilinearSimplex)

Return the physical midpoint of the 1D curvilinear element, obtained by
evaluating the mapping at `ζ = 0.5`.
"""
center(ch::CurvilinearSimplex) = pushforward(ch, 0.5)

domain(::CurvilinearSimplex{U,1,C,N,T}) where {U,C,N,T} = ReferenceSimplex{1,T,2}()

"""
    volume(ch)

Compute the physical length of the curvilinear element `ch` by integrating
the Jacobian over `[0,1]`.

A 3-point Gauss–Legendre quadrature rule is used:

∫₀¹ |dx/dζ| dζ
≈ Σ_k w_k * |dx/dζ(ζ_k)| / 2


This is the element volume used by numerical integration for BEM/FEM operators.
"""
function volume(ch::CurvilinearSimplex{U,1,C,N,T}) where {U,C,N,T}

    ξ, w = legendre(20, 0.0, 1.0)
    #ξ = (-sqrt(T(3)/T(5)), zero(T),  sqrt(T(3)/T(5)))
    #w = (T(5)/T(9),       T(8)/T(9), T(5)/T(9))
    #s = zero(T)
    #@inbounds for k in 1:3
    #    ζ = (ξ[k] + one(T)) / T(2)
    #    s += w[k] * jacobian(ch, ζ) / T(2)
    #end
    #s
    return dot(w, jacobian.(Ref(ch), ξ))
end