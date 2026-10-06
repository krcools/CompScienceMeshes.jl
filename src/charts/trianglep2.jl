# Second order triangle
struct TriangleP2{P}
    p1::P # parameter value (1,0)
    p2::P # parameter value (0,1)
    p3::P # parameter value (0,0)
    p4::P # opposite p1, parameter value (0.5, 0.5)
    p5::P # opposite p2, parameter value (0, 0.5)
    p6::P # opposite p3, parameter value (0.5, 0)
end

function coordtype(::Type{TriangleP2{P}}) where {P}
    eltype(P)
end

struct SegmentP2{P}
    p1::P # parameter value (1)
    p2::P # parameter value (0)
    p3::P # parameter value (0.5)
end

function coordtype(::Type{SegmentP2{P}}) where {P}
    eltype(P)
end

#
# Generic implementation of a Neighborhood for any chart type
#
struct Neighborhood{C,U,D,T,P}
    chart::C
    params::U
    tangents::D
    jacobian::T
    cartesian::P
end

function cartesian(p::Neighborhood) p.cartesian end

# Gmsh node ordering: p1,p2,p3 are the vertices, p4 lies on edge p1-p2,
# p5 on edge p2-p3 and p6 on edge p3-p1.
function neighborhood(ch::TriangleP2, u)
    s, t = u[1], u[2]
    l1, l2, l3 = s, t, 1 - s - t

    N = (l1*(2l1-1), l2*(2l2-1), l3*(2l3-1), 4*l1*l2, 4*l2*l3, 4*l3*l1)
    Ns = (4*l1-1, zero(l1), 1 - 4*l3, 4*l2, -4*l2, 4*(l3-l1))
    Nt = (zero(l1), 4*l2-1, 1 - 4*l3, 4*l1, 4*(l3-l2), -4*l1)

    P = (ch.p1, ch.p2, ch.p3, ch.p4, ch.p5, ch.p6)

    x  = sum(N[i]  * P[i] for i in 1:6)
    t1 = sum(Ns[i] * P[i] for i in 1:6)
    t2 = sum(Nt[i] * P[i] for i in 1:6)

    g11 = sum(t1 .* t1)
    g22 = sum(t2 .* t2)
    g12 = sum(t1 .* t2)
    j = sqrt(g11 * g22 - g12^2)

    return Neighborhood(ch, u, (t1, t2), j, x)
end


function neighborhood(ch::SegmentP2, u)
    s = u[1]
    l1, l2 = s, 1 - s

    N = (l1*(2l1-1), l2*(2l2-1), 4*l1*l2)
    Ns = (4*l1-1, 1 - 4*l2, 4*(l2-l1))

    P = (ch.p1, ch.p2, ch.p3)

    x  = sum(N[i]  * P[i] for i in 1:3)
    t1 = sum(Ns[i] * P[i] for i in 1:3)

    j = sqrt(sum(t1 .* t1))

    return Neighborhood(ch, u, (t1,), j, x)
end


function faces(ch::TriangleP2)
    return (SegmentP2(ch.p2, ch.p3, ch.p6),
            SegmentP2(ch.p3, ch.p1, ch.p4),
            SegmentP2(ch.p1, ch.p2, ch.p5))
end