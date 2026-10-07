"""
    silvester(p, x)
Compute the Silvester polynomials of order `p` at `x`.
"""
function silvester(p, x)
    T = typeof(x)

    R = zeros(T, p + 1)
    R[1] = T(1)

    for i = 0:(p-1)
        R[i + 1 + 1] = (p*x - i)/(i + 1)*R[i + 1]
    end

    return R
end

"""
    dsilvester(p, x)
Compute the derivatives of the Silvester polynomials of order `p` at `x`.
"""
function dsilvester(p, x)
    T = typeof(x)

    R = silvester(p, x)

    dR = zeros(T, p + 1)

    # (2.6)
    dR[1] = T(0)

    for i = 0:(p - 1)
        dR[i + 1 + 1] = p/(i + 1) * R[i + 1] + (p*x - i)/(i + 1) * dR[i + 1]
    end

    return dR
end

"""
    fast_silvester(p, x)
Compute the Silvester polynomials of order `p` and their derivatives at `x`.
"""
# TODO: Consider using pre-allocated memory.
function fast_silvester(p, x)
    T = typeof(x)

    R = zeros(T, p + 1)
    dR = zeros(T, p + 1)

    R[1] = T(1)
    dR[1] = T(0)

    for i = 0:(p - 1)
        R[i + 1 + 1] = (p*x - i)/(i + 1)*R[i + 1]
        dR[i + 1 + 1] = p/(i + 1) * R[i + 1] + (p*x - i)/(i + 1) * dR[i + 1]
    end

    return R, dR
end

#### Silvester polynomials up to order 10 ####
# Explicitly calculated. However, it seems that fast_silvester(p, x) 
# approach is faster despite memory allocations.
# These lines should be removed if we settle at the fast_silver approach.
@inline silv(::Val{1}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{1}, ::Val{1}, ξ::T) where {T} = ξ

@inline silv(::Val{2}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{2}, ::Val{1}, ξ::T) where {T} = 2 * ξ
@inline silv(::Val{2}, ::Val{2}, ξ::T) where {T} = 2 * ξ^2 - ξ

@inline silv(::Val{3}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{3}, ::Val{1}, ξ::T) where {T} = 3ξ
@inline silv(::Val{3}, ::Val{2}, ξ::T) where {T} = (9ξ^2 - 3ξ) / 2
@inline silv(::Val{3}, ::Val{3}, ξ::T) where {T} = (27ξ^3 - 27ξ^2 + 6ξ) / 6


@inline silv(::Val{4}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{4}, ::Val{1}, ξ::T) where {T} = 4ξ
@inline silv(::Val{4}, ::Val{2}, ξ::T) where {T} = (16ξ^2 - 4ξ) / 2
@inline silv(::Val{4}, ::Val{3}, ξ::T) where {T} = (64ξ^3 - 48ξ^2 + 8ξ) / 6
@inline silv(::Val{4}, ::Val{4}, ξ::T) where {T} = (256ξ^4 - 256ξ^3 + 96ξ^2 - 12ξ) / 24


@inline silv(::Val{5}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{5}, ::Val{1}, ξ::T) where {T} = 5ξ
@inline silv(::Val{5}, ::Val{2}, ξ::T) where {T} = (25ξ^2 - 5ξ) / 2
@inline silv(::Val{5}, ::Val{3}, ξ::T) where {T} = (125ξ^3 - 75ξ^2 + 10ξ) / 6
@inline silv(::Val{5}, ::Val{4}, ξ::T) where {T} = (625ξ^4 - 500ξ^3 + 150ξ^2 - 20ξ) / 24
@inline silv(::Val{5}, ::Val{5}, ξ::T) where {T} = (3125ξ^5 - 3125ξ^4 + 1250ξ^3 - 250ξ^2 + 20ξ) / 120


@inline silv(::Val{6}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{6}, ::Val{1}, ξ::T) where {T} = 6ξ
@inline silv(::Val{6}, ::Val{2}, ξ::T) where {T} = (36ξ^2 - 6ξ) / 2
@inline silv(::Val{6}, ::Val{3}, ξ::T) where {T} = (216ξ^3 - 108ξ^2 + 12ξ) / 6
@inline silv(::Val{6}, ::Val{4}, ξ::T) where {T} = (1296ξ^4 - 864ξ^3 + 192ξ^2 - 24ξ) / 24
@inline silv(::Val{6}, ::Val{5}, ξ::T) where {T} = (7776ξ^5 - 6480ξ^4 + 2160ξ^3 - 360ξ^2 + 30ξ) / 120
@inline silv(::Val{6}, ::Val{6}, ξ::T) where {T} = (46656ξ^6 - 46656ξ^5 + 19440ξ^4 - 4320ξ^3 + 540ξ^2 - 30ξ) / 720


@inline silv(::Val{7}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{7}, ::Val{1}, ξ::T) where {T} = 7ξ
@inline silv(::Val{7}, ::Val{2}, ξ::T) where {T} = (49ξ^2 - 7ξ) / 2
@inline silv(::Val{7}, ::Val{3}, ξ::T) where {T} = (343ξ^3 - 147ξ^2 + 14ξ) / 6
@inline silv(::Val{7}, ::Val{4}, ξ::T) where {T} = (2401ξ^4 - 1372ξ^3 + 294ξ^2 - 28ξ) / 24
@inline silv(::Val{7}, ::Val{5}, ξ::T) where {T} = (16807ξ^5 - 11760ξ^4 + 3528ξ^3 - 560ξ^2 + 35ξ) / 120
@inline silv(::Val{7}, ::Val{6}, ξ::T) where {T} = (117649ξ^6 - 100842ξ^5 + 38808ξ^4 - 8400ξ^3 + 1050ξ^2 - 42ξ) / 720
@inline silv(::Val{7}, ::Val{7}, ξ::T) where {T} = (823543ξ^7 - 823543ξ^6 + 352947ξ^5 - 88200ξ^4 + 13860ξ^3 - 1260ξ^2 + 56ξ) / 5040


@inline silv(::Val{8}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{8}, ::Val{1}, ξ::T) where {T} = 8ξ
@inline silv(::Val{8}, ::Val{2}, ξ::T) where {T} = (64ξ^2 - 8ξ) / 2
@inline silv(::Val{8}, ::Val{3}, ξ::T) where {T} = (512ξ^3 - 192ξ^2 + 16ξ) / 6
@inline silv(::Val{8}, ::Val{4}, ξ::T) where {T} = (4096ξ^4 - 2048ξ^3 + 384ξ^2 - 32ξ) / 24
@inline silv(::Val{8}, ::Val{5}, ξ::T) where {T} = (32768ξ^5 - 20480ξ^4 + 5120ξ^3 - 800ξ^2 + 48ξ) / 120
@inline silv(::Val{8}, ::Val{6}, ξ::T) where {T} = (262144ξ^6 - 196608ξ^5 + 73728ξ^4 - 16128ξ^3 + 2016ξ^2 - 72ξ) / 720
@inline silv(::Val{8}, ::Val{7}, ξ::T) where {T} = (2097152ξ^7 - 1835008ξ^6 + 802816ξ^5 - 206080ξ^4 + 32256ξ^3 - 3024ξ^2 + 84ξ) / 5040
@inline silv(::Val{8}, ::Val{8}, ξ::T) where {T} = (16777216ξ^8 - 16777216ξ^7 + 7520256ξ^6 - 2007040ξ^5 + 345600ξ^4 - 38400ξ^3 + 2688ξ^2 - 96ξ) / 40320


@inline silv(::Val{9}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{9}, ::Val{1}, ξ::T) where {T} = 9ξ
@inline silv(::Val{9}, ::Val{2}, ξ::T) where {T} = (81ξ^2 - 9ξ) / 2
@inline silv(::Val{9}, ::Val{3}, ξ::T) where {T} = (729ξ^3 - 243ξ^2 + 18ξ) / 6
@inline silv(::Val{9}, ::Val{4}, ξ::T) where {T} = (6561ξ^4 - 2916ξ^3 + 486ξ^2 - 36ξ) / 24
@inline silv(::Val{9}, ::Val{5}, ξ::T) where {T} = (59049ξ^5 - 32805ξ^4 + 8100ξ^3 - 1215ξ^2 + 54ξ) / 120
@inline silv(::Val{9}, ::Val{6}, ξ::T) where {T} = (531441ξ^6 - 354294ξ^5 + 118098ξ^4 - 24192ξ^3 + 2916ξ^2 - 90ξ) / 720
@inline silv(::Val{9}, ::Val{7}, ξ::T) where {T} = (4782969ξ^7 - 3835122ξ^6 + 1536795ξ^5 - 376200ξ^4 + 58806ξ^3 - 5292ξ^2 + 108ξ) / 5040
@inline silv(::Val{9}, ::Val{8}, ξ::T) where {T} = (43046721ξ^8 - 38742048ξ^7 + 16796160ξ^6 - 4469856ξ^5 + 777600ξ^4 - 86400ξ^3 + 6048ξ^2 - 120ξ) / 40320
@inline silv(::Val{9}, ::Val{9}, ξ::T) where {T} = (387420489ξ^9 - 387420489ξ^8 + 172186884ξ^7 - 48189030ξ^6 + 8845200ξ^5 - 1108800ξ^4 + 95040ξ^3 - 5040ξ^2 + 144ξ) / 362880

@inline silv(::Val{10}, ::Val{0}, ξ::T) where {T} = T(1)
@inline silv(::Val{10}, ::Val{1}, ξ::T) where {T} = 10ξ
@inline silv(::Val{10}, ::Val{2}, ξ::T) where {T} = (100ξ^2 - 10ξ) / 2
@inline silv(::Val{10}, ::Val{3}, ξ::T) where {T} = (1000ξ^3 - 300ξ^2 + 20ξ) / 6
@inline silv(::Val{10}, ::Val{4}, ξ::T) where {T} = (10000ξ^4 - 4000ξ^3 + 600ξ^2 - 40ξ) / 24
@inline silv(::Val{10}, ::Val{5}, ξ::T) where {T} = (100000ξ^5 - 50000ξ^4 + 10000ξ^3 - 1250ξ^2 + 50ξ) / 120
@inline silv(::Val{10}, ::Val{6}, ξ::T) where {T} = (1000000ξ^6 - 600000ξ^5 + 180000ξ^4 - 33750ξ^3 + 3750ξ^2 - 100ξ) / 720
@inline silv(::Val{10}, ::Val{7}, ξ::T) where {T} = (10000000ξ^7 - 7000000ξ^6 + 2520000ξ^5 - 588000ξ^4 + 92400ξ^3 - 8820ξ^2 + 140ξ) / 5040
@inline silv(::Val{10}, ::Val{8}, ξ::T) where {T} = (100000000ξ^8 - 80000000ξ^7 + 33600000ξ^6 - 9408000ξ^5 + 1814400ξ^4 - 241920ξ^3 + 20160ξ^2 - 160ξ) / 40320
@inline silv(::Val{10}, ::Val{9}, ξ::T) where {T} = (1000000000ξ^9 - 900000000ξ^8 + 405000000ξ^7 - 121500000ξ^6 + 27000000ξ^5 - 4320000ξ^4 + 472500ξ^3 - 31500ξ^2 + 180ξ) / 362880
@inline silv(::Val{10}, ::Val{10}, ξ::T) where {T} = (10000000000ξ^10 - 10000000000ξ^9 + 4500000000ξ^8 - 1350000000ξ^7 + 300000000ξ^6 - 50000000ξ^5 + 6000000ξ^4 - 480000ξ^3 + 24000ξ^2 - 200ξ) / 3628800

@inline dsilv(::Val{1}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{1}, ::Val{1}, ξ::T) where {T} = T(1)

@inline dsilv(::Val{2}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{2}, ::Val{1}, ξ::T) where {T} = T(2)
@inline dsilv(::Val{2}, ::Val{2}, ξ::T) where {T} = 4 * ξ - 1

@inline dsilv(::Val{3}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{3}, ::Val{1}, ξ::T) where {T} = T(3)
@inline dsilv(::Val{3}, ::Val{2}, ξ::T) where {T} = (18ξ - 3) / 2
@inline dsilv(::Val{3}, ::Val{3}, ξ::T) where {T} = (81ξ^2 - 54ξ + 6) / 6

@inline dsilv(::Val{4}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{4}, ::Val{1}, ξ::T) where {T} = T(4)
@inline dsilv(::Val{4}, ::Val{2}, ξ::T) where {T} = (32ξ - 4) / 2
@inline dsilv(::Val{4}, ::Val{3}, ξ::T) where {T} = (192ξ^2 - 96ξ + 8) / 6
@inline dsilv(::Val{4}, ::Val{4}, ξ::T) where {T} = (1024ξ^3 - 768ξ^2 + 192ξ - 12) / 24

@inline dsilv(::Val{5}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{5}, ::Val{1}, ξ::T) where {T} = T(5)
@inline dsilv(::Val{5}, ::Val{2}, ξ::T) where {T} = (50ξ - 5) / 2
@inline dsilv(::Val{5}, ::Val{3}, ξ::T) where {T} = (375ξ^2 - 150ξ + 10) / 6
@inline dsilv(::Val{5}, ::Val{4}, ξ::T) where {T} = (2500ξ^3 - 1500ξ^2 + 300ξ - 20) / 24
@inline dsilv(::Val{5}, ::Val{5}, ξ::T) where {T} = (15625ξ^4 - 12500ξ^3 + 3750ξ^2 - 500ξ + 20) / 120

@inline dsilv(::Val{6}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{6}, ::Val{1}, ξ::T) where {T} = T(6)
@inline dsilv(::Val{6}, ::Val{2}, ξ::T) where {T} = (72ξ - 6) / 2
@inline dsilv(::Val{6}, ::Val{3}, ξ::T) where {T} = (648ξ^2 - 216ξ + 12) / 6
@inline dsilv(::Val{6}, ::Val{4}, ξ::T) where {T} = (5184ξ^3 - 2592ξ^2 + 384ξ - 24) / 24
@inline dsilv(::Val{6}, ::Val{5}, ξ::T) where {T} = (38880ξ^4 - 25920ξ^3 + 6480ξ^2 - 720ξ + 30) / 120
@inline dsilv(::Val{6}, ::Val{6}, ξ::T) where {T} = (279936ξ^5 - 233280ξ^4 + 77760ξ^3 - 12960ξ^2 + 1080ξ - 30) / 720


@inline dsilv(::Val{7}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{7}, ::Val{1}, ξ::T) where {T} = T(7)
@inline dsilv(::Val{7}, ::Val{2}, ξ::T) where {T} = (98ξ - 7) / 2
@inline dsilv(::Val{7}, ::Val{3}, ξ::T) where {T} = (1029ξ^2 - 294ξ + 14) / 6
@inline dsilv(::Val{7}, ::Val{4}, ξ::T) where {T} = (9604ξ^3 - 4116ξ^2 + 588ξ - 28) / 24
@inline dsilv(::Val{7}, ::Val{5}, ξ::T) where {T} = (84035ξ^4 - 47040ξ^3 + 10584ξ^2 - 1120ξ + 35) / 120
@inline dsilv(::Val{7}, ::Val{6}, ξ::T) where {T} = (705894ξ^5 - 504210ξ^4 + 155232ξ^3 - 25200ξ^2 + 2100ξ - 42) / 720
@inline dsilv(::Val{7}, ::Val{7}, ξ::T) where {T} = (5764801ξ^6 - 4941258ξ^5 + 1764735ξ^4 - 352800ξ^3 + 41580ξ^2 - 2520ξ + 56) / 5040


@inline dsilv(::Val{8}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{8}, ::Val{1}, ξ::T) where {T} = T(8)
@inline dsilv(::Val{8}, ::Val{2}, ξ::T) where {T} = (128ξ - 8) / 2
@inline dsilv(::Val{8}, ::Val{3}, ξ::T) where {T} = (1536ξ^2 - 384ξ + 16) / 6
@inline dsilv(::Val{8}, ::Val{4}, ξ::T) where {T} = (16384ξ^3 - 6144ξ^2 + 768ξ - 32) / 24
@inline dsilv(::Val{8}, ::Val{5}, ξ::T) where {T} = (163840ξ^4 - 81920ξ^3 + 15360ξ^2 - 1600ξ + 48) / 120
@inline dsilv(::Val{8}, ::Val{6}, ξ::T) where {T} = (1572864ξ^5 - 983040ξ^4 + 294912ξ^3 - 48384ξ^2 + 4032ξ - 72) / 720
@inline dsilv(::Val{8}, ::Val{7}, ξ::T) where {T} = (14680064ξ^6 - 11010048ξ^5 + 4014080ξ^4 - 824320ξ^3 + 96768ξ^2 - 6048ξ + 84) / 5040
@inline dsilv(::Val{8}, ::Val{8}, ξ::T) where {T} = (134217728ξ^7 - 117440512ξ^6 + 45121536ξ^5 - 10035200ξ^4 + 1382400ξ^3 - 115200ξ^2 + 5376ξ - 96) / 40320


@inline dsilv(::Val{9}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{9}, ::Val{1}, ξ::T) where {T} = T(9)
@inline dsilv(::Val{9}, ::Val{2}, ξ::T) where {T} = (162ξ - 9) / 2
@inline dsilv(::Val{9}, ::Val{3}, ξ::T) where {T} = (2187ξ^2 - 486ξ + 18) / 6
@inline dsilv(::Val{9}, ::Val{4}, ξ::T) where {T} = (26244ξ^3 - 8748ξ^2 + 972ξ - 36) / 24
@inline dsilv(::Val{9}, ::Val{5}, ξ::T) where {T} = (295245ξ^4 - 131220ξ^3 + 24300ξ^2 - 2430ξ + 54) / 120
@inline dsilv(::Val{9}, ::Val{6}, ξ::T) where {T} = (3188646ξ^5 - 1771470ξ^4 + 472392ξ^3 - 72576ξ^2 + 5832ξ - 90) / 720
@inline dsilv(::Val{9}, ::Val{7}, ξ::T) where {T} = (34012203ξ^6 - 23010732ξ^5 + 7683975ξ^4 - 1504800ξ^3 + 176418ξ^2 - 10584ξ + 108) / 5040
@inline dsilv(::Val{9}, ::Val{8}, ξ::T) where {T} = (362797056ξ^7 - 310136448ξ^6 + 117573120ξ^5 - 27974400ξ^4 + 4665600ξ^3 - 518400ξ^2 + 30240ξ - 120) / 40320
@inline dsilv(::Val{9}, ::Val{9}, ξ::T) where {T} = (3874204890ξ^8 - 3486784401ξ^7 + 1205308188ξ^6 - 289134180ξ^5 + 44226000ξ^4 - 4435200ξ^3 + 285120ξ^2 - 10080ξ + 144) / 362880


@inline dsilv(::Val{10}, ::Val{0}, ξ::T) where {T} = T(0)
@inline dsilv(::Val{10}, ::Val{1}, ξ::T) where {T} = T(10)
@inline dsilv(::Val{10}, ::Val{2}, ξ::T) where {T} = (200ξ - 10) / 2
@inline dsilv(::Val{10}, ::Val{3}, ξ::T) where {T} = (3000ξ^2 - 600ξ + 20) / 6
@inline dsilv(::Val{10}, ::Val{4}, ξ::T) where {T} = (40000ξ^3 - 12000ξ^2 + 1200ξ - 40) / 24
@inline dsilv(::Val{10}, ::Val{5}, ξ::T) where {T} = (500000ξ^4 - 200000ξ^3 + 30000ξ^2 - 2500ξ + 50) / 120
@inline dsilv(::Val{10}, ::Val{6}, ξ::T) where {T} = (6000000ξ^5 - 3000000ξ^4 + 720000ξ^3 - 101250ξ^2 + 7500ξ - 100) / 720
@inline dsilv(::Val{10}, ::Val{7}, ξ::T) where {T} = (70000000ξ^6 - 42000000ξ^5 + 12600000ξ^4 - 2352000ξ^3 + 277200ξ^2 - 17640ξ + 140) / 5040
@inline dsilv(::Val{10}, ::Val{8}, ξ::T) where {T} = (800000000ξ^7 - 560000000ξ^6 + 201600000ξ^5 - 47040000ξ^4 + 7257600ξ^3 - 725760ξ^2 + 40320ξ - 160) / 40320
@inline dsilv(::Val{10}, ::Val{9}, ξ::T) where {T} = (9000000000ξ^8 - 7200000000ξ^7 + 2835000000ξ^6 - 729000000ξ^5 + 135000000ξ^4 - 17280000ξ^3 + 1417500ξ^2 - 63000ξ + 180) / 362880
@inline dsilv(::Val{10}, ::Val{10}, ξ::T) where {T} = (100000000000ξ^9 - 90000000000ξ^8 + 36000000000ξ^7 - 9450000000ξ^6 + 1800000000ξ^5 - 250000000ξ^4 + 24000000ξ^3 - 1440000ξ^2 + 48000ξ - 200) / 3628800