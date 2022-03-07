using Test

using CompScienceMeshes

for T in [Float32, Float64]
    e0 = point(T,0,0,0)
    e1 = point(T,1,0,0)
    e2 = point(T,0,1,0)
    e3 = point(T,0,0,1)

    local p = simplex(e0,e1,e2)

    @test isinside(p, (e0+e1+e2)/3) == true
    @test isinside(p, point(T, -0.1, 0.3, 0.0)) == false
    @test isinside(p, point(T, 0.3, -0.3, 0.0)) == false
    @test isinside(p, point(T, 0.6, 0.6, 0.0)) == false
    @test isinside(p, point(T, 0.4, 0.6, 0.0)) == false


    @test isinclosure(p, (e0+e1+e2)/3) == true
    @test isinclosure(p, point(T, -0.1, 0.3, 0.0)) == false
    @test isinclosure(p, point(T, 0.3, -0.3, 0.0)) == false
    @test isinclosure(p, point(T, 0.6, 0.6, 0.0)) == false
    @test isinclosure(p, point(T, 0.4, 0.6, 0.0)) == true

    s = simplex(e0,e1,e2)
    line1 = point(T,1/3,1/3,-1), point(T,1/3,1/3,+1)
    @test CompScienceMeshes.touches(s, line1) == true
    line1 = point(T,1,1,-1), point(T,1,1,+1)
    @test CompScienceMeshes.touches(s, line1) == false
end