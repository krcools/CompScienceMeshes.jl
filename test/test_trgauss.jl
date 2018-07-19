using CompScienceMeshes
using Test

ignd(x,y) = x^2 * cos(x*y)

Iref = 0.082739596439277

p1 = [0.0, 0.0]
p2 = [1.0, 0.0]
p3 = [0.0, 1.0]

N = length(CompScienceMeshes.triangleGaussW)
J = zeros(N)
for _n in 1 : N
    u, w = trgauss(_n)
    for _g in 1 : length(w)
        _p = u[1,_g]*p1 + u[2,_g]*p2 + (1-u[1,_g]-u[2,_g])*p3
        x, y = _p[1], _p[2]
        J[_n] += w[_g] * ignd(x,y)
    end
end

E = [0.3326898477547438,0.00128868370286643,0.0003505830254662883,0.0002990068871080405,0.00026769940415425264,1.4413762136992496e-7,1.8717613253691006e-7,1.8109645370453297e-12,3.4565487801639488e-12,3.0571872387445795e-12,3.3689945069086314e-12,3.3050899511418194e-12,3.3156568461898747e-12]

for _q in 1 : N
    @test abs(J[_q]-Iref)/abs(Iref) <=  E[_q]*2.0
end
