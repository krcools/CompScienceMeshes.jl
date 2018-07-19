using CompScienceMeshes
using Test

function isoriented(m::Mesh)

    @assert dimension(m) >= 0
    edges = skeleton(m, dimension(m)-1)

    D = connectivity(edges, m)
    S = (sum(D,1) .== 0)
    return all(S)
end

fn = joinpath(dirname(@__FILE__),"assets","cube2.in")
m = readmesh(fn)
@test isoriented(m)

fn = joinpath(dirname(@__FILE__),"assets","sphere2.in")
m = readmesh(fn)
@test isoriented(m)
