using Test
using LinearAlgebra
using SparseArrays
using CompScienceMeshes


for T in [Float32, Float64]
    tol = eps(T) * 10^3
    third = one(T)/3

    # test creation of a rectangle
    local m = meshrectangle(T(1.0),T(1.0),T(0.5));
    @test numvertices(m) == 9
    @test numcells(m) == 8

    # test edge detection
    m1 = skeleton(m, 1)
    @test numcells(m1) == 16

    # test volume pair construction from a list of edges
    vp = cellpairs(m,m1)
    internal_edges = count(x->x>0, vp[2,:])
    @test internal_edges == 8

    # test the relative orientation routine
    @test CompScienceMeshes.relorientation(index(1,2,3),index(1,2,3,4)) == +4
    @test CompScienceMeshes.relorientation(index(2,3,4),index(1,2,3,4)) == -1
    @test CompScienceMeshes.relorientation(index(1,3,2),index(1,2,3,4)) == -4
    @test CompScienceMeshes.relorientation(index(2,4,3),index(1,2,3,4)) == +1

    o, e1, e2, e3 = euclidianbasis(3)
    cell = simplex(o, e1, e2)

    A  = volume(cell)
    @test abs(A - T(0.5)) < tol

    mp = neighborhood(cell,[third, third]);
    r = cartesian(mp)
    @test norm(r - point(third, third, zero(T))) < tol

    # repeat the test but now on a random cell
    randpoint() = point(2*rand(T)-1, 2*rand(T)-1, 2*rand(T)-1)
    v1 = randpoint()
    v2 = randpoint()
    v3 = randpoint()
    cell = simplex(v1, v2, v3)

    mp = neighborhood(cell,[third, third]);
    r = cartesian(mp)

    @test norm(third*(v1 + v2 + v3) - r) < tol
    @test abs(dot(cross(r-v1, r-v2), r-v3)) < tol

    # test that the normal is directed outwards
    m = meshcircle(T(1.0), T(2π/51))
    #for i in 1:numcells(m)
    for cl in cells(m)
        #p = chart(m, cells(m,i))
        v = vertices(m, cl)
        p = simplex(v)
        _c = (p.vertices[1] + p.vertices[2]) / 2
        @test dot(_c, p.normals[1]) > 0
    end
end

## test cellpairs on non-oriented meshes
using CompScienceMeshes
using Test

for T in [Float32, Float64]
    p1 = point(T,1,0,0)
    p2 = point(T,0,1,0)
    p3 = point(T,-1,0,0)
    p4 = point(T,0,-1,0)

    i1 = index(1,2,3)
    i2 = index(1,4,2)

    local m = Mesh([p1,p2,p3,p4], [i1,i2])
    e = skeleton(m,1)

    @test numcells(m) == 2
    @test numcells(e) == 5

    cps = cellpairs(m,e)

    #Test unstructed mesh rectangle
    r1 = meshrectangle(T(1.0),T(1.0),T(1/1))
    # r2 = meshrectangle(1.0,1.0,1/1,structured=false)
    r2 = readmesh(joinpath(@__DIR__,"assets/rectangle1.in"),T=T)

    @test numcells(r1) == 2
    @test numcells(r2) == 4
    @test numvertices(r2) == 5
    @test numcells(boundary(r2)) == 4
    @test numcells(CompScienceMeshes.interior(r2)) == 4
    @test numcells(boundary(r1)) == numcells(boundary(r2))

    #Test mesh cuboid parts
    deltaT = T(1/1)
    box = meshcuboid(T(1.0),T(1.0),T(1.0),deltaT)
    top = meshcuboid(T(1.0),T(1.0),T(1.0),deltaT,physical="TopPlate")
    sides = meshcuboid(T(1.0),T(1.0),T(1.0),deltaT,physical="SidePlates")
    bot = meshcuboid(T(1.0),T(1.0),T(1.0),deltaT,physical="BottomPlate")
    obox = meshcuboid(T(1.0),T(1.0),T(1.0),deltaT,physical="OpenBox")

    @test numcells(boundary(box)) == 0
    @test numcells(boundary(top)) == numcells(boundary(bot))  == 4
    @test numvertices(box) == 14
    @test numcells(box) == 24
    @test numcells(skeleton(box,1)) == 36
    @test numcells(top) + numcells(bot) + numcells(sides) == numcells(box)
    @test numcells(obox) == numcells(box) - numcells(top)

    e1 = skeleton(top,1)
    E1 = skeleton(box,1)

    e1_cells = cells(e1)
    E1_cells = cells(E1)

    # cols = vec([6 16 17 13 18 19 20 21])
    # Z11 = sparse(collect(1:8),cols,ones(8),numcells(e1),numcells(E1))
    if T == Float64
        S11 = CompScienceMeshes.embedding(e1,E1)
        @test size(S11) == (length(e1), length(E1))
        vals = nonzeros(S11)
        rows = rowvals(S11)
        for col in axes(S11,2)
            for k in nzrange(S11,col)
                row = rows[k]
                val = vals[k]
                val == 0 && continue
                @test e1_cells[row] == E1_cells[col]
            end
        end
    end
# @test S11 == Z11
end