using Base.Test
using CompScienceMeshes

sp = CompScienceMeshes.SphereChart(point(1,2,3), 1.0)
for L in 0:10
    @eval @test sum(w for (p,w) in quadpoints(sp, $L)) ≈ 4π
end

function Ylm(K, r)

    @assert K >= 2 "Multipole order needs to be equal or larger than 2."

    r = r / norm(r)
    x = r[1]; y = r[2]; z = r[3]

    T = eltype(r)
    Y = zeros(Complex{T}, K*K)

    Y[1] = 1/√4π

    facphip =  x + im*y;
    facphim = -x + im*y;
    for m ∈ 1:(K-1)
        Y[(m+1)*(m+1)] = -Y[m*m] * facphip * sqrt((m+0.5)/m)
        Y[m*m+1] = -Y[(m-1)*(m-1)+1] * facphim * sqrt((m+0.5)/m)
    end

    Y[3] = √3 * z * Y[1]
    for l ∈ 2:(K-1)
        Y[(l+1)*(l+1)-1] = sqrt(2l+1) * z * Y[l*l]
        Y[l*l+2] = sqrt(2*l+1) * z * Y[(l-1)*(l-1)+1]
    end

    for l ∈ 2:(K-1)
        for m ∈ (2-l):(l-2)
            Y[l*(l+1)+m+1] =
                sqrt((4*l^2-1)/(l^2-m^2)) * z * Y[l*(l-1)+m+1] -
                sqrt(((2*l+1)*((l-1)^2-m^2))/((l^2-m^2)*(2*l-3))) * Y[(l-1)*(l-2)+m+1]
        end
    end

    Y
end


K = 3
I = zeros(K^2)
for (p,w) in quadpoints(sp, 20)
    Y = Ylm(K, cartesian(p))
    for l in 0:K-1
        for m in -l:l
            i = l*(l+1)+m+1
            I[i] += w * real(Y[i] * conj(Y[i]))
        end
    end
end

@test all(isapprox.(I, 1.0))
