using CompScienceMeshes
using StaticArrays
using Test

for T in [Float32, Float64]
    local p1 = point(T, 1.0, 0.0, 0.0)
    local p2 = point(T, 10.0, 0.0, 1.0)
    local p3 = point(T, 0.0, 2.0, 5.0)

    ch = simplex(p1,p2,p3)

    dom = domain(ch)
    q, u = quadpoints(dom, 7)
    @which neighborhood(ch, q[1])

    pw = CompScienceMeshes.quadpoints(ch, 7);

    #t1 = sum(w)
    t1 = sum(q[2] for q in pw)
    t2 = T(0.5)*norm(cross(p1-p3, p2-p3))
    @test 2 * (t1-t2) / (t1+t2) + 1  ≈ 1

    p1 = point(T, 1.0, 0.0, 0.0)
    p2 = point(T, 2.0, 0.0, 0.0)
    seg = simplex(p1,p2)

    wpvs = quadpoints(x->T(cartesian(x)[1]^2), [seg], (20,))
    global r = T(0.0)
    for wpv in wpvs[1,1]
        w = wpv.weight
        v = wpv.value
        r += w * v
    end

    @test r ≈ T(1/3*(2^3-1^3))
end
