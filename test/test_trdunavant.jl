using CompScienceMeshes
using Test
using StaticArrays

#= The integrand used to test the quadrature is given by

f(x,y) = x^o + y^o + x^(floor(o/2)) * y^(ceil(o/2)),

where o is the order of the polynomial is used for testing.

Each quadrature rule 'r' is supposed to integrate
the corresponding polynomial resulting by 'o:=r', exactly.

Thus, the functions used are

functions = [
    "x¹  + y¹  +  y",
    "x²  + y²  + xy",
    "x³  + y³  + xy²",
    "x⁴  + y⁴  + x²y²",
    "x⁵  + y⁵  + x²y³",
    "x⁶  + y⁶  + x³y³",
    "x⁷  + y⁷  + x³y⁴",
    "x⁸  + y⁸  + x⁴y⁴",
    "x⁹  + y⁹  + x⁴y⁵",
    "x¹⁰ + y¹⁰ + x⁵y⁵",
    "x¹¹ + y¹¹ + x⁵y⁶",
    "x¹² + y¹² + x⁶y⁶",
    "x¹³ + y¹³ + x⁶y⁷",
    "x¹⁴ + y¹⁴ + x⁷y⁷",
    "x¹⁵ + y¹⁵ + x⁷y⁸",
    "x¹⁶ + y¹⁶ + x⁸y⁸",
    "x¹⁷ + y¹⁷ + x⁸y⁹",
    "x¹⁸ + y¹⁸ + x⁹y⁹",
    "x¹⁹ + y¹⁹ + x⁹y¹⁰",
    "x²⁰ + y²⁰ + x¹⁰y¹⁰",
]
=#

@testitem "accuracy" begin
using StaticArrays
ignd(x, y, order::Int) = (x^order + y^order + x^(ceil(order/2)) * y^(floor(order/2)))

N = length(CompScienceMeshes.trianglequadDunavantW)

for T in @SVector [Float32, Float64]

    local p1 = T.([0.0, 0.0])
    local p2 = T.([1.0, 0.0])
    local p3 = T.([0.0, 1.0])

    # Analytical Results
    local ana = @SVector (T[
        0.5,                    # 1       /   2,
        0.2083333333333333400,  # 5       /   24,
        0.1166666666666666700,  # 7       /   60,
        0.0722222222222222200,  # 13      /   180,
        0.05,                   # 1       /   20,
        0.0366071428571428600,  # 41      /   1120,
        0.0281746031746031750,  # 71      /   2520,
        0.0223809523809523800,  # 47      /   2100,
        0.0182539682539682550,  # 23      /   1260,
        0.0151815776815776810,  # 505     /   33264,
        0.0128343878343878340,  # 925     /   72072,
        0.0109949574235288520,  # 1849    /   168168,
        0.0095265845265845270,  # 3433    /   360360,
        0.0083345473970473980,  # 1373    /   164736,
        0.0073535125005713240,  # 12871   /   1750320,
        0.0065362016342408500,  # 25741   /   3938220,
        0.0058480734951323185,  # 8314191 /   1421697420,
        0.0052632120201779640,  # 97241   /   18475600,
        0.0047619305359243440,  # 3879897 /   814773960,
        0.0043290160444677760,  # 123171  /   28452424,
    ])

    J = zeros(N)
    for _n in 1 : N
        u, w = trgauss(TriangleQuadDunavant(_n))
        for _g in eachindex(w)
            _p = u[1,_g]*p1 + u[2,_g]*p2 + (1-u[1,_g]-u[2,_g])*p3
            x, y = _p[1], _p[2]
            J[_n] += w[_g] * ignd(x, y, _n)
        end
    end

    # Tolerance factor of 200 for Float64 necessary
    T == Float32 ? tol = T(1) : tol = T(200)

    for _q in 1 : N
        #@show abs(J[_q]-ana[_q])/abs(ana[_q])
        @test abs(J[_q]-ana[_q])/abs(ana[_q]) <= tol*eps(T)
    end
end
end