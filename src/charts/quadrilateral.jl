struct Quadrilateral{P}
    p1::P
    p2::P
    p3::P
    p4::P
end

function coordtype(q::Quadrilateral{P}) where {P} eltype(P) end
function vertices(q::Quadrilateral) SVector(q.p1,q.p2,q.p3,q.p4) end

function cartesian(quad::Quadrilateral, u)
    return quad.p1 + u[1] * (quad.p2 - quad.p1) + u[2] * (quad.p4 - quad.p1) + u[1] * u[2] * (quad.p3 - quad.p4 + quad.p1 - quad.p2)
end

function tangents(quad::Quadrilateral, u)
    aux = quad.p3 - quad.p4 + quad.p1 - quad.p2
    ∂ru = quad.p2 - quad.p1 + u[2] * aux
    ∂rv = quad.p4 - quad.p1 + u[1] * aux
    return hcat(∂ru, ∂rv)
end

function normal(quad::Quadrilateral, u)
    Du = tangents(quad, u)
    normalize(Du[:,1] × Du[:,2])
end

function jacobian(quad::Quadrilateral, u)
    Du = tangents(quad, u)
    g = Du' * Du
    √det(g)
end

function faces(ch::Quadrilateral)
    return SVector(
        simplex(ch.p1,ch.p2),
        simplex(ch.p2,ch.p3),
        simplex(ch.p3,ch.p4),
        simplex(ch.p4,ch.p1),)
end


struct Neighborhood{C,P,Q,T,J,N}
    chart::C
    parametric::P
    cartesian::Q
    tangents::T
    jacobian::J
    normal::N
end


function neighborhood(quad::Quadrilateral, u)
    c = cartesian(quad, u)
    t = tangents(quad, u)
    q = t[:,1] × t[:,2]
    j = norm(q)
    n = normalize(q)
    Neighborhood(quad, u, c, t, j, n)
end


function parametric(p::Neighborhood) p.parametric end
function cartesian(p::Neighborhood) p.cartesian end
function tangents(p::Neighborhood) p.tangents end
function tangents(p::Neighborhood, i::Int) p.tangents[:,i] end
function normal(p::Neighborhood) p.normal end
function jacobian(p::Neighborhood) p.jacobian end


function neighborhood_lazy(quad::Quadrilateral, u) neighborhood(quad, u) end


@testitem "quadrilateral chart" begin
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)
    p = neighborhood(quad, point(0.5, 0.5))
    @show cartesian(p)
    @test cartesian(p) ≈ point(0.5,0.5,0.0)
    @test tangents(p, 1) ≈ point(1,0,0)
    @test tangents(p, 2) ≈ point(0,1,0)
end





struct RefQuadrilateral{T} end
domain(quad::Quadrilateral{P}) where {P} = RefQuadrilateral{eltype(P)}()
neighborhood(quad::RefQuadrilateral, u) = SVector(u)
neighborhood(ch::RefQuadrilateral, u::AbstractVector) = SVector{length(u)}(u)

function faces(ch::RefQuadrilateral)
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    return SVector(
        simplex(p1,p2),
        simplex(p2,p3),
        simplex(p3,p4),
        simplex(p4,p1),)
end


function quadpoints(ch::RefQuadrilateral{T}, rule) where {T}

    U1, W1 = legendre(rule, zero(T), one(T))
    U2, W2 = legendre(rule, zero(T), one(T))

    [(neighborhood(ch, (u1,u2)), w1*w2) for (u1,w1) in zip(U1,W1) for (u2,w2) in zip(U2,W2)]
end

@testitem "RefQuadrilateral: quadpoints" begin
    refchart = CompScienceMeshes.RefQuadrilateral{Float64}()
    qps = quadpoints(refchart, 5)
    fn = p -> 1.0
    I = sum(w*fn(p) for (p,w) in qps)
    @test I ≈ 1.0
end

@testitem "Quadrilateral: quadpoints" begin
    p1 = point(0,0,0)
    p2 = point(2,0,0)
    p3 = point(2,3,0)
    p4 = point(0,3,0)
    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)
    qps = quadpoints(quad, 5)
    fn = p -> 1.0
    I = sum(w*fn(p) for (p,w) in qps)
    @test I ≈ 6.0
end

function permute_vertices(q::Quadrilateral, I)
    verts = vertices(q)[I]
    return Quadrilateral(verts...)
end

function center(q::Quadrilateral)
    T = coordtype(q)
    h = T(0.5)
    return neighborhood(q, (h,h))
end