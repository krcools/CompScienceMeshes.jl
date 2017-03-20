immutable WeightPointValue{W,P,V}
  weight::W
  point::P
  value::V
end


function quadpoints(chart::ReferenceSimplex{1}, rule)
    u, w = legendre(rule, 0.0, 1.0)
    [neighborhood(chart, u[:,i]) for i in eachindex(w)], w
end

function quadpoints(chart::ReferenceSimplex{2}, rule)
    u, w = trgauss(rule)
    [neighborhood(chart, u[:,i]) for i in eachindex(w)], w
end


function quadpoints(chart, rule)
    P, V = quadpoints(domain(chart), rule)
    Q = [neighborhood(chart,p) for p in P]
    W = [jacobian(q)*v for (q,v) in zip(Q,V)]
    return Q, W
end


"""
    quadpoints(refspace, charts, rules)

Computed a matrix of vectors containing (weight, point, value) triples that can
be used in numerical integration over the elements described by the charts. Internally,
this method used `quadpoints(chart, rule)` to retrieve the points and weights for
a certain quadrature rule over `chart`.
"""
function quadpoints(f, charts, rules)

    P = pointtype(eltype(charts))
    p, w = quadpoints(charts[1], rules[1])
    W = eltype(w)
    V = typeof(f(p[1]))

    WPV = WeightPointValue{W,P,V}

    qd = Array{Vector{WPV}}(length(rules), length(charts))
    for j in eachindex(charts)
        for i in eachindex(rules)
            p, w = quadpoints(charts[j], rules[i])
            qd[i,j] = Vector{WPV}(length(w))
            for k in eachindex(w)
                #wk = w[k] * jacobian(p[k])
                qd[i,j][k] = WeightPointValue(w[k],p[k],f(p[k]))
            end
        end
    end

  qd
end
