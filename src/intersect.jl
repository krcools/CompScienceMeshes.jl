export intersection

"""
    function intersect{U,C,T}(p1::FlatCellNM{U,1,C,2,T}, p2::FlatCellNM{U,1,C,2,T})

Compute the intersection of two overlapping segments.
"""
function intersection{U,C,T}(p1::FlatCellNM{U,1,C,2,T}, p2::FlatCellNM{U,1,C,2,T})

    tol = sqrt(eps(T))
    PT = eltype(p1.vertices)

    W = [u for u in p2.vertices]
    clipconvex!(W, p1.vertices[1],  p1.tangents[1])
    clipconvex!(W, p1.vertices[2], -p1.tangents[1])

    return patch(W, Val{1})
end

function clipconvex!(W, v, m)
    m2 = dot(m,m)
    for i in 1:length(W)
        s = dot(W[i]-v,m)
        s = min(0, s)
        W[i] = v + s / m2 * m
    end
end

"""
  function intersection{U,C,T}(p1::FlatCellNM{U,2,C,3,T}, p2::FlatCellNM{U,2,C,3,T})

Returns the intersection of triangles p1 and p2. ATTENTION: The implementation
as is assumes that either p1 is a subset of p2 or the other way around. This
case is sufficient to compute matrices for integral operators on combinations
of RT and BC spaces defined subordinate to the same mesh.
"""
function intersection{U,C,T}(p1::FlatCellNM{U,2,C,3,T}, p2::FlatCellNM{U,2,C,3,T})
  b = true
  for v in p2.vertices
    if !isinclosure(p1,v)
      b = false
      break
    end
  end

  b == true ? p2 : p1
end
