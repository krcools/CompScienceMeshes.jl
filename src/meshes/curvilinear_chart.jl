#### Implementation background, todos etc
# We implement the standard Lagrange interpolated curvilinear elements
# For the arbitrary order code, we use the Silvester polynomials as described in
# [1] R. D. Graglia and A. F. Peterson, Higher-order techniques in computational electromagnetics.
# TODO:
# 1) Extend to 3D simplices.
# 2) Add as many specializations as needed
# 3) Rethink what MeshPointNM should contain as information

"""
    CurvilinearSimplex{U,D,C,N,T}

A curvilinear simplex (currently specialized to 1D) embedded in `U`-dimensional
space, represented by `N` interpolation vertices and `N` parametric nodes.

# Fields
- `vertices :: SVector{N, SVector{U,T}}`  
    Coordinates of the interpolation/control points in physical space.

The type parameters are:
- `U`: ambient (universe) dimension,
- `D`: topological dimension,
- `C`: codimension (= `U-D`),
- `N`: number of local nodes,
- `T`: coordinate type (e.g. `Float64`).
- `O`: element order.

This is the geometric building block used for curvilinear 1D meshes.
"""
struct CurvilinearSimplex{U,D,C,N,T,O} <: AbstractSimplex{U,D}
    vertices::SVector{N,SVector{U,T}}end

"""
    nodes(s::CurvilinearSimplex{U,1})

Return the physical coordinates of the two endpoints of a 1D curvilinear element,
identified by the parametric values `ζ = 1` (left) and `ζ = 0` (right) under the
barycentric convention `ζ = λ_left`.

# Returns
`(x_left, x_right) :: (SVector{U,T}, SVector{U,T})`
"""
function nodes(s::CurvilinearSimplex{U,1}) where {U}
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

barytocart(ch::CurvilinearSimplex{U,1}, u::Real) where {U} = barytocart(ch, SVector{1}(u))
barytocart(ch::CurvilinearSimplex{U,2}, u::Real, v::Real) where {U} = barytocart(ch, SVector{2}(u, v))
barytocart(ch::CurvilinearSimplex{U,3}, u::Real, v::Real, w::Real) where {U} = barytocart(ch, SVector{3}(u, v, w))

cartesian(ch::AbstractSimplex{U,1}, u::SVector{1,<:Real})  where {U} = barytocart(ch, u)
parametric(ch::AbstractSimplex{U,1}, u::SVector{1,<:Real}) where {U} = u

function barytocart(ch::CurvilinearSimplex{U,1,C,N,T,O}, bary::AbstractVector) where {U,C,N,T,O}
    u = bary[1]

    s = zero(SVector{U,T})
    @inbounds for (r, (i, j)) in enumerate(gmsh_line_index_to_tuple[O])
        s += ch.vertices[r] * silvester(O, u)[i + 1] * silvester(O, 1 - u)[j + 1]
    end

    return s
end

# barytocart specializations for 1D curvilinear elements using Silvester interpolation
# Orders 1 through 10, expressions simplified and verified such that u=1 yields v1

function barytocart(ch::CurvilinearSimplex{U,1,C,N,T,1}, bary::AbstractVector) where {U,C,N,T}
    u = bary[1]
    return ch.vertices[1] * u + ch.vertices[2] * (1 - u)
end

function barytocart(ch::CurvilinearSimplex{U,1,C,N,T,2}, bary::AbstractVector) where {U,C,N,T}
    u = bary[1]
    u2 = u*u

    return ch.vertices[1] * (2u2 - u) +
           ch.vertices[2] * (2u2 - 3u + 1) +
           ch.vertices[3] * (-4u2 + 4u)
end

function barytocart(ch::CurvilinearSimplex{U,2,C,N,T,O}, bary) where {U,C,N,T,O}
    u = bary[1]
    v = bary[2]

    s = zero(SVector{U,T})
    @inbounds for (r, (i, j, k)) in enumerate(gmsh_triangle_index_to_tuple[O])
        s += ch.vertices[r] * 
            silvester(O, u)[i + 1] *
            silvester(O, v)[j + 1] *
            silvester(O, 1 - u - v)[k + 1]
    end

    return s
end

function barytocart(ch::CurvilinearSimplex{U,2,C,N,T,1}, bary) where {U,C,N,T}
    u, v = bary
    w = 1 - u - v

    return ch.vertices[1] * u +
           ch.vertices[2] * v +
           ch.vertices[3] * w
end

#=
function barytocart(ch::CurvilinearSimplex{U,2,C,N,T,2}, bary) where {U,C,N,T}
    u, v = bary
    w = 1 - u - v

    u2 = u * u
    v2 = v * v
    w2 = w * w

    return ch.vertices[1] * (2u2 - u) +
           ch.vertices[2] * (2v2 - v) +
           ch.vertices[3] * (2w2 - w) +
           ch.vertices[4] * (4u*v) +
           ch.vertices[5] * (4v*w) +
           ch.vertices[6] * (4u*w)
end

function barytocart(ch::CurvilinearSimplex{U,2,C,N,T,3}, bary) where {U,C,N,T}
    u, v = bary
    w = 1 - u - v

    u2 = u * u
    u3 = u2 * u
    v2 = v * v
    v3 = v2 * v
    w2 = w * w
    w3 = w2 * w

    return ch.vertices[1] * (10u3 - 15u2 + 6u) +
           ch.vertices[2] * (10v3 - 15v2 + 6v) +
           ch.vertices[3] * (10w3 - 15w2 + 6w) +
           ch.vertices[4] * (60u2*v - 60u*v) +
           ch.vertices[5] * (60u*v2 - 60u*v) +
           ch.vertices[6] * (60v2*w - 60v*w) +
           ch.vertices[7] * (60v*w2 - 60v*w) +
           ch.vertices[8] * (60u*w2 - 60u*w) +
           ch.vertices[9] * (60u2*w - 60u*w) +
           ch.vertices[10] * (60u*v*w)
end
=#

"""
    unitarybase(ch::CurvilinearSimplex, bary)

Compute the unitary base vectors as defined in (3.38) in Graglia & Peterson.
"""
function unitarybase(ch::CurvilinearSimplex{U,1,C,N,T,O}, bary) where {U,C,N,T,O}
    u = bary[1]

    s = zero(SVector{U,T})

    Ru, dRu = fast_silvester(O, u)
    Rv, dRv = fast_silvester(O, 1 - u)

    @inbounds for (r, (i, j)) in enumerate(gmsh_line_index_to_tuple[O])

        # Adapted (2.44) in Graglia & Peterson with ξ₁ = u and ξ₂ = 1 - u
        s += ch.vertices[r] * (dRu[i + 1] * Rv[j + 1] - Ru[i + 1] * dRv[j + 1])
    end

    return SVector{1}(s)
end

function unitarybase(ch::CurvilinearSimplex{U,1,C,N,T,1}, bary) where {U,C,N,T}
    u = bary[1]

    # Adapted (2.44) in Graglia & Peterson with ξ₁ = u and ξ₂ = 1 - u
    s = ch.vertices[1] - ch.vertices[2]

    return SVector{1}(s)
end

function unitarybase2(ch::CurvilinearSimplex{U,1,C,N,T,2}, bary) where {U,C,N,T}
    u = bary[1]

    # Adapted (2.44) in Graglia & Peterson with ξ₁ = u and ξ₂ = 1 - u
    s = ch.vertices[1] * (+4u - 1) +
        ch.vertices[2] * (+4u - 3) +
        ch.vertices[3] * (-8u + 4)

    return SVector{1}(s)
end

# TODO: Speed-up gains possible by specialization to explicit polynomials
# do at least for lowest orders to avoid memory allocations
function unitarybase(ch::CurvilinearSimplex{U,2,C,N,T,O}, bary) where {U,C,N,T,O}
    u = bary[1]
    v = bary[2]

    Ru, dRu = fast_silvester(O, u)
    Rv, dRv = fast_silvester(O, v)
    Rw, dRw = fast_silvester(O, 1 - u - v)

    𝓵¹ = zero(SVector{U,T})
    𝓵² = zero(SVector{U,T})
    @inbounds for (r, (i, j, k)) in enumerate(gmsh_triangle_index_to_tuple[O])

        𝓵¹ += ch.vertices[r] * (
                dRu[i + 1] * Rv[j + 1] * Rw[k + 1] -
                Ru[i + 1] * Rv[j + 1] * dRw[k + 1]
            )

        𝓵² += ch.vertices[r] * (
                Ru[i + 1] * dRv[j + 1] * Rw[k + 1] -
                Ru[i + 1] * Rv[j + 1] * dRw[k + 1]
            )
    end
    
    return SVector(𝓵¹, 𝓵²)
end

# TODO: Trying to avoid allocations. Even for lowest order however this is much slower. Why???
function unitarybase_experimental(ch::CurvilinearSimplex{U,2,C,N,T,1}, bary) where {U,C,N,T}
    u = bary[1]
    v = bary[2]

    Ru = @SVector[silv(Val(1), Val(i), u) for i=0:1]
    Rv = @SVector[silv(Val(1), Val(i), v) for i=0:1]
    Rw = @SVector[silv(Val(1), Val(i), 1 - u - v) for i=0:1]

    dRu = @SVector[dsilv(Val(1), Val(i), u) for i=0:1]
    dRv = @SVector[dsilv(Val(1), Val(i), v) for i=0:1]
    dRw = @SVector[dsilv(Val(1), Val(i), 1 - u - v) for i=0:1]

    𝓵¹ = zero(SVector{U,T})
    𝓵² = zero(SVector{U,T})
    @inbounds for (r, (i, j, k)) in enumerate(gmsh_triangle_index_to_tuple[1])

        𝓵¹ += ch.vertices[r] * (
                dRu[i + 1] * Rv[j+1] * Rw[k+1] -
                Ru[i + 1] * Rv[j+1] * dRw[k+1]
            )

        𝓵² += ch.vertices[r] * (
                Ru[i + 1] * dRv[j+1] * Rw[k+1] -
                Ru[i + 1] * Rv[j+1] * dRw[k+1]
            )
    end
    
    return SVector(𝓵¹, 𝓵²)
end


function unitarybase2(ch::CurvilinearSimplex{U,2,C,N,T,7}, bary) where {U,C,N,T}
    u = bary[1]
    v = bary[2]

    Ru = @SVector[silv(Val(7), Val(i), u) for i=0:7]
    Rv = @SVector[silv(Val(7), Val(i), v) for i=0:7]
    Rw = @SVector[silv(Val(7), Val(i), 1 - u - v) for i=0:7]

    dRu = @SVector[dsilv(Val(7), Val(i), u) for i=0:7]
    dRv = @SVector[dsilv(Val(7), Val(i), v) for i=0:7]
    dRw = @SVector[dsilv(Val(7), Val(i), 1 - u - v) for i=0:7]

    𝓵¹ = zero(SVector{U,T})
    𝓵² = zero(SVector{U,T})
    @inbounds for (r, (i, j, k)) in enumerate(gmsh_triangle_index_to_tuple[7])

        𝓵¹ += ch.vertices[r] * (
                dRu[i + 1] * Rv[j+1] * Rw[k+1] -
                Ru[i + 1] * Rv[j+1] * dRw[k+1]
            )

        𝓵² += ch.vertices[r] * (
                Ru[i + 1] * dRv[j+1] * Rw[k+1] -
                Ru[i + 1] * Rv[j+1] * dRw[k+1]
            )
    end
    
    return SVector(𝓵¹, 𝓵²)
end

function chart(
    m::CompScienceMeshes.CurvilinearMesh{U,D1,N,T,O},
    faceid::Int
) where {U,D1,N,T,O}

    vindices = m.faces[faceid]
    X = m.vertices[vindices]

    CompScienceMeshes.CurvilinearSimplex{U,D1-1,U-1,N,T,O}(X)
end

function normal(ch::CurvilinearSimplex{U,2}, bary) where {U}
    t = unitarybase(ch, bary)

    nt = cross(t[1], t[2])

    return nt ./ norm(nt)
end


"""
    jacobian(ch, ζ)

Return the scalar Jacobian determinant.
"""
jacobian(ch::CurvilinearSimplex{U,1}, u) where {U} = norm(unitarybase(ch, u))

function jacobian(ch::CurvilinearSimplex{U,2}, bary) where {U} 
    t = unitarybase(ch, bary)

    return norm(cross(t[1], t[2]))
end

"""
    center(ch::CurvilinearSimplex)

Return the physical midpoint of the 1D curvilinear element, obtained by
evaluating the mapping at `ζ = 0.5`.
"""
center(ch::CurvilinearSimplex{U,1}) where {U} = barytocart(ch, SVector(0.5))
center(ch::CurvilinearSimplex{U,2}) where {U} = barytocart(ch, SVector(0.5, 0.5))

domain(::CurvilinearSimplex{U,1,C,N,T}) where {U,C,N,T} = ReferenceSimplex{1,T,2}()
domain(::CurvilinearSimplex{U,2,C,N,T}) where {U,C,N,T} = ReferenceSimplex{2,T,3}()


"""
    volume(ch)

Compute the physical length of the curvilinear element `ch` by integrating
the Jacobian over `[0,1]`.
"""
function volume(ch::CurvilinearSimplex{U,1,C,N,T}) where {U,C,N,T}

    ξ, w = legendre(20, 0.0, 1.0)

    return dot(w, jacobian.(Ref(ch), ξ))
end

"""
    volume(ch)

Compute the physical length of the curvilinear element `ch` by integrating
the Jacobian over `[0,1]`.
"""
# TODO: volume should be saved in the simplex type. However, Jacobian requires the existence...
function volume(ch::CurvilinearSimplex{U,2,C,N,T}) where {U,C,N,T}

    ξ, w = legendre(20, 0.0, 1.0)

    vol = T(0)

    for i in eachindex(ξ)
        for j in eachindex(ξ)
            vol += w[i] * w[j] * jacobian(ch, SVector(ξ[i], ξ[j]*(1 - ξ[i]))) * (1 - ξ[i])
        end
    end

    return vol
end