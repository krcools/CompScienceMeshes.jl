#export Simplex


# U: the dimension of the universe
# D: the dimension of the manifold
# N: the number of vertices
# T: the type of the coordinates
# C: the complimentary dimension (should always be U-D)
struct Simplex{U,D,C,N,T}
    vertices::SVector{N,SVector{U,T}}
    tangents::SVector{D,SVector{U,T}}
    normals::SVector{C,SVector{U,T}}
    volume::T
end
normal(t::Simplex{3,2,1,3,<:Number}) = t.normals[1]
normal(t::Simplex{3,2,1,3,<:Number}, u) = t.normals[1] 
dimtype(splx::Simplex{U,D}) where {U,D} = Val{D}
"""
    permute_simplex(simplex,permutation)

Permutation is a Vector v which sets the v[i]-th vertex at the i-th place.

Return Simplex with permuted vertices list, tangents are recalculated, normal is kept the same
"""
function permute_vertices(s::Simplex{3,2,1,3,T},permutation::Union{Vector{P},SVector{P}}) where {T,P}
    vert = SVector{3,SVector{3,T}}(s.vertices[permutation])
    simp = simplex(vert)
    Simplex(simp.vertices,simp.tangents,s.normals,simp.volume)
end

"""
    flip_normal(simplex)

Flips the normal of the simplex. Only on triangles embedded in 3D space
"""
flip_normal(t::Simplex{3,2,1,3,<:Number}) = Simplex(t.vertices,t.tangents,-t.normals,t.volume)

"""
    flip_normal(simplex, sign)

Flips the normal of the simplex if sign is -1. Only on triangles embedded in 3D space
"""
function flip_normal(t::Simplex{3,2,1,3,<:Number},sign::Int)
    sign == 1 && return t
    return flip_normal(t)
    end
export flip_normal

# tangents(s::Simplex{3,2,1,3,<:Number},i) = s.tangents[i]
"""
    coordtype(simplex)

Return coordinate type used by simplex.
"""
coordtype(::Type{Simplex{U,D,C,N,T}}) where {U,D,C,N,T} = T
coordtype(p::Simplex) = coordtype(typeof(p))


"""
    volume(simplex)

Return the volume of the simplex.
"""
volume(p::Simplex) = p.volume


"""
A tuple of points, aka an interval behaves trivially like a chart
"""
volume(x::Tuple{T,T}) where {T} = norm(x[2]-x[1])


"""
    dimension(simplex)

Return the manifold dimension of the simplex.
"""
dimension(::Type{Simplex{U,D,C,N,T}}) where {U,D,C,N,T} = D

"""
    dimension(simplex)

Return the manifold dimension of the simplex.
"""
dimension(p::Simplex) = dimension(typeof(p))


"""
    length(simplex)

Returns the number of vertices (equals dimension + 1)
"""
Base.length(p::Simplex) = dimension(typeof(p))+1


"""
  universedimension(p)

Return the dimension of the universe in which `p` is embedded.
"""
universedimension(::Type{Simplex{U,D,C,N,T}}) where {U,D,C,N,T} = U
universedimension(p::Simplex) = universedimension(typeof(p))


"""
    getindex(simplex, I)
    simplex[I]

Get the vertices at index I (scalar or array) defining the simplex
"""
getindex(p::Simplex, I::Union{Number,SVector,Array}) = p.vertices[I]



"""
    simplex(vertices)
    simplex(v1, v2, ...)
    simplex(vertices, Val{D})

Build a D-dimensional simplex. The vertices can be passed in
an array (static or dynamic), or supplied separately. If the
length of the array is not part of its type, the speed of the
construction can be improved by supplying an extra Val{D}
argument. In case it is not clear from the context whether
the vertex array is dynamically or statically sized, use the
third form as it will not incur notable performance hits.

Note that D is the dimension of the simplex, i.e. the number
of vertices supplied minus one.
"""
@generated function simplex(vertices::SVector{D1,P}) where {D1,P}
    U = length(P)
    D = D1 - 1
    C = U-D
    T = eltype(P)
    xp1 =:(())
    for i in 1:D
        push!(xp1.args, :(vertices[$i]-vertices[end]))
    end
    xp2 = :(SVector{$D,P}($xp1))
    quote
        tangents = $xp2
        normals, volume = _normals(tangents, Val{$C})
        Simplex(vertices, tangents, normals, $T(volume))
    end
end

simplex(vertices::SVector{N,T}...) where {N,T<:Number} = simplex(SVector((vertices...,)))

@generated function simplex(vertices, ::Type{Val{D}}) where D
    P = eltype(vertices)
    D1 = D + 1
    xp = :(())
    for i in 1:D1
        push!(xp.args, :(vertices[$i]))
    end
    :(simplex(SVector{$D1,$P}($xp)))
end

# normal(s::Simplex{3,2,1,3,T}) where {T} = normalize(cross(s[1]-s[3],s[2]-s[3]))


function _normals(tangents, ::Type{Val{1}})
    PT = eltype(tangents)
    D  = length(tangents)
    T  = eltype(PT)

    n = zeros(T,D+1)
    b = Array{T}(undef,D,D)

    for i in 1:D+1
        fill!(b, zero(T))
        for j in 1:D
            for k in 1:i-1
                b[k,j] = tangents[j][k]
            end
            for k in i:D
                b[k,j] = tangents[j][k+1]
            end
        end
        n[i] = (-1)^(i-1) * det(b)
    end

    n *= (-1)^D / norm(n)
    normals = SVector{1,PT}(PT(n))

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) /  factorial(D)

    return normals, volume

end

function _normals(tangents::SVector{2,SVector{3,T}}, ::Type{Val{1}}) where {T}

    n = tangents[1] × tangents[2]
    l = norm(n)

    P = SVector{3,T}
    SVector{1,P}(n/l), 0.5*l
end

function _normals(tangents::SVector{2,SVector{2,T}}, ::Type{Val{0}}) where {T}

    t = tangents[1]
    s = tangents[2]
    v = (t[1]*s[2] - t[2]*s[1])/2
    # n[3] = tangents[1] × tangents[2]
    # l = norm(n)

    P = SVector{2,T}
    SVector{0,P}(), v
end

function _normals(tangents, ::Type{Val{C}}) where C
    PT = eltype(tangents)
    D  = length(tangents)
    U = length(PT)
    T  = eltype(PT)

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) / factorial(D)

    # Fix this. This function needs to become gneerated
    normals = SVector{C,PT}(PT[zero(PT) for i in 1:C])

    return normals, volume
end


function _normals(tangents::SVector{1,SVector{3,T}} where {T}, ::Type{Val{2}})
    P = eltype(tangents)
    normals = SVector{2,P}(zero(P), zero(P))
    # volume = dot(cross(tangents[1], tangents[2]), tangents[3]) / 6
    volume = norm(tangents[1])
    return normals, volume
end

function _normals(tangents::SVector{3,SVector{3,T}} where {T}, ::Type{Val{0}})
    P = eltype(tangents)
    normals = SVector{0,P}()
    volume = dot(cross(tangents[1], tangents[2]), tangents[3]) / 6
    return normals, volume
end



"""
    barytocart(simplex, uv)

Returns the point in the simplex with barycentric coordinates uv
"""
function barytocart(mani::Simplex, u)
    r = last(mani.vertices)
    for i in 1 : dimension(mani)
        ti = mani.tangents[i]
        ui = u[i]
        #r += mani.tangents[i] * u[i]
        r += ti * ui
    end
    return r
end

function cartesian(ch::Simplex, u) barytocart(ch, u) end


"""
    carttobary(simplex, point) -> barycoords

Compute the barycentric coordinates on 'simplex' of 'point'.
"""
function carttobary(p::Simplex{U,D,C,N,T}, cart) where {U,D,C,N,T}

    G = [dot(p.tangents[i], p.tangents[j]) for i in 1:D, j in 1:D]
    #w = [dot(p.tangents[i], cart - p.vertices[end]) for i in 1:D]

    o = p.vertices[end]
    w = [sum(t[j]*(cart[j]-o[j]) for j in 1:length(cart)) for t in p.tangents]

    u = G \ w

    return SVector{D}(u)
end

const _combs24 = [SVector{2}(c) for c in combinations(@SVector[1,2,3,4],2)]
const _edgeidx24 = [relorientation(c, SVector(1,2,3,4)) for c in _combs24]
const _combs24_pos = [(p > 0 ? _combs24[i] : reverse(_combs24[i])) for (i,p) in enumerate(_edgeidx24)]
const _edgeidx24_abs = abs.(_edgeidx24)

function edges(s::Simplex{3,3})
    T = eltype(eltype(s.vertices))
    P = Simplex{3,1,2,2,T}
    Edges = Vector{P}(undef, length(_combs24_pos))
    for (i,c) in zip(_edgeidx24_abs, _combs24_pos)
        edge = simplex(s.vertices[c[1]], s.vertices[c[2]])
        Edges[i] = edge
    end
    return Edges
end


function edges(s::Simplex{3,2})
    return [
        simplex(s.vertices[2], s.vertices[3]),
        simplex(s.vertices[3], s.vertices[1]),
        simplex(s.vertices[1], s.vertices[2])]
end

function faces(c::Simplex{3,2})
    @SVector[
        simplex(c[2], c[3]),
        simplex(c[3], c[1]),
        simplex(c[1], c[2]),
    ]
end

function faces(c::CompScienceMeshes.Simplex{3,3,0,4,T}) where {T}
    @SVector[
        simplex(c[4],c[3],c[2]),
        simplex(c[1],c[3],c[4]),
        simplex(c[1],c[4],c[2]),
        simplex(c[1],c[2],c[3])
    ]
end

function faces(ch::Simplex{3,2}, ::Type{Val{0}})
    simplex.(ch.vertices)
end

function subcharts(c::Simplex{U,2}, ::Type{Val{1}}) where {U}
    T = coordtype(c)
    d = (
        point(T, 1, 0),
        point(T, 0, 1),
        point(T, 0, 0),
    )
    tp = (
        (simplex(c[2], c[3]), simplex(d[2], d[3])),
        (simplex(c[3], c[1]), simplex(d[3], d[1])),
        (simplex(c[1], c[2]), simplex(d[1], d[2])),
    )
    return tp
end

@testitem "subcharts" begin

    ch = simplex(
        point(2,0,-1),
        point(0,3,-1),
        point(0,0,-1))

    fc = CompScienceMeshes.subcharts(ch, Val{1})
    @test length(fc) == 3
    @test volume(fc[1][1]) ≈ 3
    @test volume(fc[2][1]) ≈ 2
    @test volume(fc[3][1]) ≈ sqrt(13)

end


"""
    ReferenceSimplex{Dimension, CoordType, NumVertices}

This domain is defined to bootstrap the quadrature generation strategy. The generic definition
of numquads on a chart pulls back to the domain. For a limit set of reference domains, explicit
quadrature rules are defined. The weights and points are then pushed forward to the configuaration
space element over which integration is desired.

For more details see the implementation in quadpoints.jl
"""
struct ReferenceSimplex{D,T,N}
    # N: number of defining points
    # D: dimension of the simplex
    # T: type of the coordinates
    # simplex::Simplex{D,D,0,N,T}
end

# function ReferenceSimplex{D,T,N}() where {D,T,N}
#     P = SVector{D,T}[]
#     for i in 1:D
#         a = zeros(T,D)
#         a[i] = 1
#         p = SVector{D,T}(a)
#         push!(P,p)
#     end

#     o = zero(SVector{D,T})
#     push!(P,o)
#     ReferenceSimplex{D,T,N}(simplex(P...))
# end

function faces(c::ReferenceSimplex{2})
    p = vertices(c)
    return SVector(
        simplex(p[2], p[3]),
        simplex(p[3], p[1]),
        simplex(p[1], p[2]))
end

function faces(ch::ReferenceSimplex{2}, ::Type{Val{0}})
    simplex.(vertices(ch))
end

barytocart(ch::ReferenceSimplex, u) = u # barytocart(ch.simplex, u)
carttobary(ch::ReferenceSimplex, p) = SVector{2}(p[1:2]) # carttobary(ch.simplex, p)

function domain(ch::Simplex{U,D,C,N,T}) where {U,D,C,T,N} ReferenceSimplex{D,T,N}() end
function domain(ch::ReferenceSimplex) ch end

neighborhood(ch::ReferenceSimplex, u) = SVector(u)
neighborhood(ch::ReferenceSimplex, u::AbstractVector) = SVector{length(u)}(u)

function vertices(ch::ReferenceSimplex{1,T}) where {T}
    SVector(
        point(T,1),
        point(T,0),
    )
end

function vertices(ch::ReferenceSimplex{2,T}) where {T}
    SVector(
        point(T,1,0),
        point(T,0,1),
        point(T,0,0))
end

function permute_vertices(ch::ReferenceSimplex, I)
    V = vertices(ch)
    return simplex(SVector(V[I[1]], V[I[2]], V[I[3]]))
end

# points act as neighborhoods for the identity chart
parametric(u) = u
cartesian(u) = u
jacobian(u) = one(eltype(u))
function tangents(u::SVector{N}) where {N}
    SMatrix{N,N}(LinearAlgebra.I)
end
tangents(u::SVector{N}, i::Int) where {N} = tangents(u)[:,i]


# """
#     tangents(splx)

# Returns a matrix whose columns are the tangents of the simplex `splx`.
# """
# tangents(splx::Simplex) = hcat((splx.tangents)...)
tangents(splx::Simplex, u) = hcat((splx.tangents)...)

# vertices(splx::Simplex) = hcat((splx.vertices)...)
function vertices(s::Simplex) s.vertices end

"""
    verticeslist(simplex)

Returns the vertices as a list.
"""
verticeslist(splx::Simplex) = splx.vertices
export verticeslist