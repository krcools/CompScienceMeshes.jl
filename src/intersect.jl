"""
    intersect(chartA, chartB) -> [chart1, chart2, ...]

Compute the intersection of two charts of equal dimension.

Compute an array of charts such that the disjoint union of those charts produces the intersection of the two charts provided as inputs. In particular the sum of integrals over the returned charts will equal the integral over the intersection of the two given charts.
"""
function intersection(p1::Simplex{U,1,C,2,T}, p2::Simplex{U,1,C,2,T}) where {U,C,T}

    tol = sqrt(eps(T))
    PT = eltype(p1.vertices)

    W = [u for u in p2.vertices]
    clipconvex!(W, p1.vertices[1],  p1.tangents[1])
    clipconvex!(W, p1.vertices[2], -p1.tangents[1])

    # For consistency an array needs to be return. In higher
    # dimensions the intersection could be the union of
    # multiple simplices.
    [simplex(W[1], W[2])]
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
    intersection(triangleA, triangle B)

The output inherits the orientation from the first operand
"""
function intersection(p1::Simplex{U,2,C,3}, p2::Simplex{U,2,C,3}) where {U,C}
  pq = sutherlandhodgman(p1.vertices, p2.vertices)
  return [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
end

function intersection(p1::Simplex{3,2,1,3}, p2::Simplex{3,2,1,3};tol=eps())
    pq = sutherlandhodgman(p1.vertices, p2.vertices)
    nonoriented_simplexes = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_simplexes = nonoriented_simplexes[volume.(nonoriented_simplexes).>tol]
    signs = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(normal(p1)))))
    flip_normal.(nonoriented_simplexes,signs)
  end


function intersection(p1::Simplex{U,3,C,4}, p2::Simplex{U,3,C,4}) where {U,C}
  @assert overlap(p1,p2)
  volume(p1) <= volume(p2) ? [p1] : [p2]
end

function intersection_keep_clippings(p1::Simplex{3,2,1,3}, p2::Simplex{3,2,1,3};tol=1e3 * eps(coordtype(p1)))
    pq, cls = sutherlandhodgman_keep_clippings(p1.vertices, p2.vertices)
    # @show size(cls)
    nonoriented_simplexes = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_simplexes = nonoriented_simplexes[volume.(nonoriented_simplexes).>tol]

    nonoriented_clippings = map(cls) do cl
            [ simplex(cl[1], cl[i], cl[i+1]) for i in 2:length(cl)-1]
    end
    nonoriented_clippings = map(nonoriented_clippings) do cl
        cl[volume.(cl) .> tol]
    end

    n1 = normal(p1)
    signs = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(n1))))
    oriented_simplices = flip_normal.(nonoriented_simplexes,signs)
    oriented_clippings = map(nonoriented_clippings) do cl
        signs = Int.(sign.(dot.(normal.(cl), Ref(n1))))
        flip_normal.(cl, signs)
    end
    return oriented_simplices, oriented_clippings
end

"""
    intersection2(triangleA, triangleB)

returns intersection in which the first operand inherits orientation from first argument
and second argument inherits orientation of second argument.
"""
function intersection2(p1::Simplex, p2::Simplex)
    a = intersection(p1,p2)
    return [(i,i) for i in a]
end

function intersection2(p1::Simplex{3,2,1,3}, p2::Simplex{3,2,1,3}; tol=eps())
    pq = sutherlandhodgman(p1.vertices, p2.vertices)
    nonoriented_simplexes = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_simplexes = nonoriented_simplexes[volume.(nonoriented_simplexes).>tol]
    signs1 = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(normal(p1)))))
    
    signs2 = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(normal(p2)))))
    return [(flip_normal(s,signs1[i]),flip_normal(s,signs2[i])) for (i,s) in enumerate(nonoriented_simplexes)]
end
export intersection2


"""
    intersectline(a,b,p,q)

Computes the intersection of the lines (in a 2D space) defined
by points [a,b] and [p,q]
"""
function intersectlines(a,b,p,q)

    P = typeof(a)

    d1 = det(@SMatrix [a[1] a[2]; b[1] b[2]])
    d2 = det(@SMatrix [p[1] p[2]; q[1] q[2]])

    id = one(eltype(a))
    d1x = det(@SMatrix [a[1] id; b[1] id])
    d1y = det(@SMatrix [a[2] id; b[2] id])

    d2x = det(@SMatrix [p[1] id; q[1] id])
    d2y = det(@SMatrix [p[2] id; q[2] id])

    den = det(@SMatrix [d1x d1y; d2x d2y])
    @assert !isinf(1/den)

    P(
        det(@SMatrix [d1 d1x; d2 d2x]) / den,
        det(@SMatrix [d1 d1y; d2 d2y]) / den,
    )

end



"""
    inside(p,a,b)

Determines is p is on the interior side of the segment b of the boundary,
assuming that the boundary is oriented counter-clockwise.
"""
function leftof(p,a,b)

    tol = sqrt(eps(eltype(a)))
    ap = @SVector [ p[1]-a[1], p[2]-a[2] ]
    ab = @SVector [ b[1]-a[1], b[2]-a[2] ]
    ap[1]*ab[2] - ap[2]*ab[1] <= tol ? true : false

end

"""
    sutherlandhodgman2d(subject,clipper)

Computes the intersection of the coplanar triangles
subject and clipper.
"""
function sutherlandhodgman2d(subject,clipper)

    PT = eltype(subject)

    clipped = Array(copy(subject))
    sizehint!(clipped, 8)

    input = Array(copy(clipped))
    sizehint!(input, 8)

    b = last(clipper)
    for a in clipper

        resize!(input, length(clipped))
        copyto!(input, clipped)
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

function sutherlandhodgman2d_keep_clippings(subject,clipper)

    clipped = Vector(copy(subject))
    clippings = Vector{typeof(clipped)}()

    b = last(clipper)
    for a in clipper

        input = deepcopy(clipped)
        resize!(clipped,0)
        clipping = similar(clipped, 0)

        if !isempty(input)
            q = last(input)
            for p in input
                if leftof(p, b, a)
                    if !leftof(q, b, a)
                        # q is right, p is left
                        i = intersectlines(b,a,q,p)
                        push!(clipped, i)
                        push!(clipping, i)
                        push!(clipped, p)
                    else
                        # q is left, p is left
                        push!(clipped, p)
                    end

                else # p is right
                    if leftof(q, b, a)
                        # q is left, p is right
                        i = intersectlines(b,a,q,p)
                        push!(clipped, i)

                        push!(clipping, i)
                        push!(clipping, p)
                    else
                        push!(clipping, p)
                    end
                end

                q = p
            end
        end
        push!(clippings, clipping)
        b = a
    end

    return clipped, clippings
end


"""
    sutherlandhodgman(subject, clipper)

Compute the intersection of two coplanar triangles, potentially
embedded in a higher dimensional space.
"""
function sutherlandhodgman(subject, clipper)
    triangle = simplex(clipper, Val{2})
    subject2d = [carttobary(triangle,p) for p in subject]
    for p in subject2d for x in p @assert !isinf(x) end end
    clipper2d = [carttobary(triangle,q) for q in clipper]
    for p in clipper2d for x in p @assert !isinf(x) end end
    clipped2d = sutherlandhodgman2d(subject2d, clipper2d)
    clipped = [barytocart(triangle,q) for q in clipped2d ]
end

"""
    sutherlandhodgman_keep_clippings(subject, clipper)

Compute both the intersection and the setminus of two convex polygons `subject` and
`clipper`. The polygons should be supplied as vectors of points. The intersection is
returned as a vector of points and the clippings as a vector of vectors of points. 
"""
function sutherlandhodgman_keep_clippings(subject, clipper)
    triangle = simplex(clipper, Val{2})
    subject2d = [carttobary(triangle,p) for p in subject]
    clipper2d = [carttobary(triangle,q) for q in clipper]
    clipped2d, clippings2d = sutherlandhodgman2d_keep_clippings(subject2d, clipper2d)
    clipped = [barytocart(triangle,q) for q in clipped2d ]
    clippings = [ [ barytocart(triangle,q) for q in clipping2d] for clipping2d in clippings2d ]
    return clipped, clippings
end
