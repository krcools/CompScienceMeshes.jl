include("TriangleDunavant.jl")
include("TriangleGauss.jl")

export trgauss
export TriangleQuadDunavant
export TriangleQuadLegacy

abstract type TriangleQuadrature end

struct TriangleQuadDunavant{I} <: TriangleQuadrature
    order::I
end

struct TriangleQuadLegacy{I} <: TriangleQuadrature
    order::I
end

"""
    trgauss(n) -> (u,w)

Returns the n-th triangle quadrature rule. Returns a Matrix u of size (Q,2)
with Q the number of quadrature points and a Vector w of size (Q,) containing
the quadrature weights.
"""
function trgauss(n)
    return trgauss(TriangleQuadLegacy(n))
end

function trgauss(trqd::TriangleQuadDunavant)
    # n = trqd.order
    # @assert 1 <= n <= length(trianglequadDunavantA)
    # u = copy(transpose([trianglequadDunavantA[n] trianglequadDunavantB[n]]))
    # w = trianglequadDunavantW[n]
    # return u, w/2

    n = trqd.order
    @assert 1 <= n <= length(trianglequadDunavantA)
    A = trianglequadDunavantA[n]
    B = trianglequadDunavantB[n]
    T = eltype(A)
    AB = [A, B]
    u = T[AB[i][j] for i in (1,2), j in 1:length(A)]
    w = trianglequadDunavantW[n]
    return u, w/2
end

function trgauss(trqd::TriangleQuadLegacy)
    n = trqd.order
    @assert 1 <= n <= length(trianglequadGaussA)
    A = trianglequadGaussA[n]
    B = trianglequadGaussB[n]
    T = eltype(A)
    AB = [A, B]
    u = T[AB[i][j] for i in (1,2), j in 1:length(A)]
    w = trianglequadGaussW[n]
    return u, w/2
end