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

#function intersectlines(p,q)
function intersectlines(a,b,p,q)

    d1 = det(@fsa [a[1] a[2]; b[1] b[2]])
    d2 = det(@fsa [p[1] p[2]; q[1] q[2]])

    id = one(eltype(a))
    d1x = det(@fsa [a[1] id; b[1] id])
    d1y = det(@fsa [a[2] id; b[2] id])

    d2x = det(@fsa [p[1] id; q[1] id])
    d2y = det(@fsa [p[2] id; q[2] id])

    den = det(@fsa [d1x d1y; d2x d2y])

    Point(
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

    ap = @fsa [ p[1]-a[1], p[2]-a[2] ]
    ab = @fsa [ b[1]-a[1], b[2]-a[2] ]
    ap[1]*ab[2] - ap[2]*ab[1] <= 0 ? true : false

end

export sutherlandhodgman
function sutherlandhodgman(subject,clipper)

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

function sutherlandhodgman3d(subject, clipper)
    T = eltype(eltype(subject))

    subject2d = zeros(Point{2,T},length(subject))
    for i in 1:length(subject)
        subject2d[i] = carttobary(subject[i])
    end

    clipper2d = zeros(Point{2,T},length(clipper))
    for i in 1:length(clipper)
        clipper2d[i] = carttobary(clipper[i])
    end

    isct = surtherlandhodgman(subject2d, clipper2d)

end
