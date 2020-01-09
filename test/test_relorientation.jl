using CompScienceMeshes
using Test

using StaticArrays

tet = @SVector[:a,:b,:c,:d]

@test CompScienceMeshes.relorientation(@SVector[:a,:b,:c], tet) == 4
@test CompScienceMeshes.relorientation(@SVector[:c,:a,:b], tet) == 4
@test CompScienceMeshes.relorientation(@SVector[:b,:c,:a], tet) == 4
@test CompScienceMeshes.relorientation(@SVector[:a,:c,:b], tet) == -4
@test CompScienceMeshes.relorientation(@SVector[:b,:a,:c], tet) == -4
@test CompScienceMeshes.relorientation(@SVector[:c,:b,:a], tet) == -4

@test CompScienceMeshes.relorientation(@SVector[:d,:c,:b], tet) == 1
@test CompScienceMeshes.relorientation(@SVector[:b,:d,:c], tet) == 1
@test CompScienceMeshes.relorientation(@SVector[:c,:b,:d], tet) == 1
@test CompScienceMeshes.relorientation(@SVector[:c,:d,:b], tet) == -1
@test CompScienceMeshes.relorientation(@SVector[:b,:c,:d], tet) == -1
@test CompScienceMeshes.relorientation(@SVector[:d,:b,:c], tet) == -1

@test CompScienceMeshes.relorientation(@SVector[:a,:c,:d], tet) == 2
@test CompScienceMeshes.relorientation(@SVector[:d,:a,:c], tet) == 2
@test CompScienceMeshes.relorientation(@SVector[:c,:d,:a], tet) == 2
@test CompScienceMeshes.relorientation(@SVector[:c,:a,:d], tet) == -2
@test CompScienceMeshes.relorientation(@SVector[:d,:c,:a], tet) == -2
@test CompScienceMeshes.relorientation(@SVector[:a,:d,:c], tet) == -2

@test CompScienceMeshes.relorientation(@SVector[:a,:d,:b], tet) == 3
@test CompScienceMeshes.relorientation(@SVector[:b,:a,:d], tet) == 3
@test CompScienceMeshes.relorientation(@SVector[:d,:b,:a], tet) == 3
@test CompScienceMeshes.relorientation(@SVector[:d,:a,:b], tet) == -3
@test CompScienceMeshes.relorientation(@SVector[:b,:d,:a], tet) == -3
@test CompScienceMeshes.relorientation(@SVector[:a,:b,:d], tet) == -3
