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

    # For consistency an array needs to be return. In higher
    # dimensions the intersection could be the union of
    # multiple simplices.
    return [patch(W, Val{1})]
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
  function intersection{U,C}(p1::FlatCellNM{U,2,C,3}, p2::FlatCellNM{U,2,C,3})

Returns the intersection of triangles p1 and p2. ATTENTION: The implementation
as is assumes that either p1 is a subset of p2 or the other way around. This
case is sufficient to compute matrices for integral operators on combinations
of RT and BC spaces defined subordinate to the same mesh.
"""
function intersection{U,C}(p1::FlatCellNM{U,2,C,3}, p2::FlatCellNM{U,2,C,3})
  # b = true
  # for v in p2.vertices
  #   if !isinclosure(p1,v)
  #     b = false
  #     break
  #   end
  # end
  #
  # b == true ? p2 : p1

  pq = sutherlandhodgman(p1.vertices, p2.vertices)
  return [ patch([pq[1], pq[i], pq[i+1]], Val{2}) for i in 2:length(pq)-1 ]
end



"""
    intersectline(a,b,p,q)

Computes the intersection of the lines (in a 2D space) defined
by points [a,b] and [p,q]
"""
function intersectlines(a,b,p,q)

    P = typeof(a)

    d1 = det(@fsa [a[1] a[2]; b[1] b[2]])
    d2 = det(@fsa [p[1] p[2]; q[1] q[2]])

    id = one(eltype(a))
    d1x = det(@fsa [a[1] id; b[1] id])
    d1y = det(@fsa [a[2] id; b[2] id])

    d2x = det(@fsa [p[1] id; q[1] id])
    d2y = det(@fsa [p[2] id; q[2] id])

    den = det(@fsa [d1x d1y; d2x d2y])

    @assert !isinf(1/den)
    #
    # if isnan(d2)
    #     @show p
    #     @show q
    # end
    # @assert !isnan(d1)
    # @assert !isnan(d2)
    # @assert !isnan(d1x)
    # @assert !isnan(d2x)
    # @assert !isnan(d1y)
    # @assert !isnan(d2y)
    # @assert !isnan(den)

    P(
        det(@fsa [d1 d1x; d2 d2x]) / den,
        det(@fsa [d1 d1y; d2 d2y]) / den,
    )

end



"""
    inside(p,a,b)

Determines is p is on the interior side of the segment b of the boundary,
assuming that the boundary is oriented counter-clockwise.
"""
function leftof(p,a,b)

    tol = sqrt(eps(eltype(a)))
    ap = @fsa [ p[1]-a[1], p[2]-a[2] ]
    ab = @fsa [ b[1]-a[1], b[2]-a[2] ]
    ap[1]*ab[2] - ap[2]*ab[1] <= tol ? true : false

end

export sutherlandhodgman2d

"""
    sutherlandhodgman2d(subject,clipper)

Computes the intersection of the coplanar triangles
subject and clipper.
"""
function sutherlandhodgman2d(subject,clipper)

    PT = eltype(subject)

    clipped = copy(subject)
    sizehint!(clipped, 8)

    input = copy(clipped)
    sizehint!(input, 8)

    b = last(clipper)
    for a in clipper

        resize!(input, length(clipped))
        copy!(input, clipped)
        resize!(clipped, 0)

        isempty(input) || (q = last(input))

        for p in input
            if leftof(p, b, a)
                if !leftof(q, b, a)
                    ist = intersectlines(b,a,q,p)
                    push!(clipped, ist)
                end

                push!(clipped, p)

            elseif leftof(q, b, a)

                ist = intersectlines(b,a,q,p)
                push!(clipped, ist)
            end

            q = p
        end
        b = a
    end

    return clipped
end

export sutherlandhodgman

"""
    sutherlandhodgman(subject, clipper)

Compute the intersection of two coplanar triangles, potentially
embedded in a higher dimensional space.
"""
function sutherlandhodgman(subject, clipper)

    triangle = patch(clipper, Val{2})
    subject2d = [carttobary(triangle,p) for p in subject]
    for p in subject2d for x in p @assert !isinf(x) end end
    clipper2d = [carttobary(triangle,q) for q in clipper]
    for p in clipper2d for x in p @assert !isinf(x) end end
    clipped2d = sutherlandhodgman2d(subject2d, clipper2d)
    clipped = [barytocart(triangle,q) for q in clipped2d ]

end
