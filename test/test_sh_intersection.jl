using FixedSizeArrays
using CompScienceMeshes
using Base.Test

p = @fsa [0.0, 1.0]
q = @fsa [2.0, 1.0]
r = @fsa [2.0, 2.0]
a = @fsa [1.0, 0.0]
b = @fsa [2.0, 0.0]
c = @fsa [1.0, 2.0]

@test CompScienceMeshes.leftof(p, a,b) == true
@test CompScienceMeshes.intersectlines(p,r, a,c) == [1.0, 1.5]

α = sutherlandhodgman([p,q,r], [a,b,c])
@test length(α) == 4
@test α[1] == [1.2, 1.6]
@test α[2] == [1.0, 1.5]
@test α[3] == [1.0, 1.0]
@test α[4] == [1.5, 1.0]

@time sutherlandhodgman([p,q,r], [a,b,c])
#@time for i in 1:450_000 sutherlandhodgman([p,q,r], [a,b,c]) end
