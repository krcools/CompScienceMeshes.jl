"""
    point(xs...)

Create point of default type (double precision coordinates)
"""
point(xs...) = point(Float64, xs...)

"""
    point(type, xs...)

Create point of default type and supplied precision for the coordinates
"""
@generated function point(T::Type,xs...)
    D = length(xs)
    xp = :(SVector{$D,T}())
    for d in 1:D
        push!(xp.args, :(xs[$d]))
    end
    return xp
end


"""
    index(i1, i2, ...) -> ids

Create a tuple of vertex indices.
"""
index(is...) = SVector{length(is),Int}(is...)



"""
  euclidian_basis(type, dim)

Returns the origin and default unit vectors for Euclidian space of dimension dim
"""
function euclidianbasis(T::DataType, dim)
  P = SVector{dim,T}
  id = eye(dim)
  r = P[ P(id[:,i]...) for i in 1:dim ]
  z = P(zeros(T,dim)...)
  return unshift!(r, z)
end
