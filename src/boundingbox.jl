"""
Returns the boundingbox of a patch in terms of its center and halfsize.

    function boundingbox{U,D,C,N,T}(p::Simplex{U,D,C,N,T}) -> center, halfsize
"""
function boundingbox(p::Simplex{U,D,C,N,T}) where {U,D,C,N,T}

    # ll = minimum(p.vertices)
    # ur = maximum(p.vertices)

    ll = first(p.vertices); for v in p.vertices; ll = min.(v, ll); end
    ur = first(p.vertices); for v in p.vertices; ll = max.(v, ll); end

    center = (ll + ur) / 2
    halfsize = maximum(ur - center)

    return center, halfsize
end

function boundingbox(V)

    # ll = minimum(V)
    # ur = maximum(V)

    ll = first(V); for v ∈ V; ll = min.(v, ll); end
    ur = first(V); for v ∈ V; ur = max.(v, ur); end

    c = (ll + ur) / 2
    s = maximum(ur - c)

    return c, s
end


boundingbox(p::SVector{N,T}) where {N,T<:Number} = p, zero(eltype(p))
