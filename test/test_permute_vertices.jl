using CompScienceMeshes
using Test

splx = simplex(
    point(0,0,0),
    point(1,0,0),
    point(0,1,0)
);

@test normal(splx) ≈ point(0,0,1)

splx1 = CompScienceMeshes.permute_vertices(splx, [2,1,3])
@test normal(splx1) ≈ point(0,0,1)

splx2 = CompScienceMeshes.flip_normal(splx)
@test normal(splx2) ≈ point(0,0,-1)

for i in 1:2
    @test splx2.tangents[i] ≈ splx.tangents[i]
end

splx3 = simplex(
    point(1,0.5,0),
    point(-1,0.5,0),
    point(-1,2.5,0)
)

isct = CompScienceMeshes.intersection2(splx, splx3)
@test length(isct) == 1

splx4, splx5 = isct[1]
@test volume(splx4) ≈ volume(splx5)
@test normal(splx4) ≈ -normal(splx5)

