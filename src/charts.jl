#export Simplex


# U: the dimension of the universe
# D: the dimension of the manifold
# N: the number of vertices
# T: the type of the coordinates
# C: the complimentary dimension (should always be U-D)
immutable Simplex{U,D,C,N,T}
    vertices::SVector{N,SVector{U,T}}
    tangents::SVector{D,SVector{U,T}}
    normals::SVector{C,SVector{U,T}}
    volume::T
end


"""
    coordtype(simplex)

Return coordinate type used by simplex.
"""
coordtype{U,D,C,N,T}(::Type{Simplex{U,D,C,N,T}}) = T
coordtype(p::Simplex) = coordtype(typeof(p))


"""
    volume(simplex)

Return the volume of the simplex.
"""
volume(p::Simplex) = p.volume


"""
A tuple of points, aka an interval behaves trivially like a chart
"""
volume{T}(x::Tuple{T,T}) = norm(x[2]-x[1])


"""
    dimension(simplex)

Return the manifold dimension of the simplex.
"""
dimension{U,D,C,N,T}(::Type{Simplex{U,D,C,N,T}}) = D
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
universedimension{U,D,C,N,T}(::Type{Simplex{U,D,C,N,T}}) = U
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
@generated function simplex{D1,P}(vertices::SVector{D1,P})
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
        Simplex(vertices, tangents, normals, volume)
    end
end

simplex(vertices...) = simplex(SVector((vertices...)))

@generated function simplex{D}(vertices, ::Type{Val{D}})
    P = eltype(vertices)
    D1 = D + 1
    xp = :(())
    for i in 1:D1
        push!(xp.args, :(vertices[$i]))
    end
    :(simplex(SVector{$D1,$P}($xp)))
end


function _normals(tangents, ::Type{Val{1}})
    PT = eltype(tangents)
    D  = length(tangents)
    T  = eltype(PT)

    n = zeros(T,D+1)
    b = Array{T}(D,D)

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
    normals = SVector{1,PT}([PT(n)])

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) /  D

    return normals, volume

end




function _normals{C}(tangents, ::Type{Val{C}})
    PT = eltype(tangents)
    D  = length(tangents)
    U = length(PT)
    T  = eltype(PT)

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) / D

    # Fix this. This function needs to become gneerated
    normals = SVector{C,PT}([zero(PT) for i in 1:C])

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


"""
    carttobary(simplex, point) -> barycoords

Compute the barycentric coordinates on 'simplex' of 'point'.
"""
function carttobary{U,D,C,N,T}(p::Simplex{U,D,C,N,T}, cart)

    G = [dot(p.tangents[i], p.tangents[j]) for i in 1:D, j in 1:D]
    #w = [dot(p.tangents[i], cart - p.vertices[end]) for i in 1:D]

    o = p.vertices[end]
    w = [sum(t[j]*(cart[j]-o[j]) for j in 1:length(cart)) for t in p.tangents]

    u = G \ w

    return SVector{D}(u)
end


"""
    ReferenceSimplex{Dimension, CoordType, NumVertices}

This domain is defined to bootstrap the quadrature generation strategy. The generic definition
of numquads on a chart pulls back to the domain. For a limit set of reference domains, explicit
quadrature rules are defined. The weights and points are then pushed forward to the configuaration
space element over which integration is desired.

For more details see the implementation in quadpoints.jl
"""
immutable ReferenceSimplex{D,T,N}
    # N: number of defining points
    # D: dimension of the simplex
    # T: type of the coordinates
    simplex::Simplex{D,D,0,N,T}
end

function (::Type{ReferenceSimplex{D,T,N}}){D,T,N}()
    P = SVector{D,T}[]
    for i in 1:D
        a = zeros(T,D)
        a[i] = 1
        p = SVector{D,T}(a)
        push!(P,p)
    end

    o = zero(SVector{D,T})
    push!(P,o)

    ReferenceSimplex(simplex(P...))
end


barytocart(ch::ReferenceSimplex, u) = barytocart(ch.simplex, u)
carttobary(ch::ReferenceSimplex, p) = carttobary(ch.simplex, p)

domain{U,D,C,T,N}(ch::Simplex{U,D,C,N,T}) = ReferenceSimplex{D,T,N}()
neighborhood(ch::ReferenceSimplex, u) = u
