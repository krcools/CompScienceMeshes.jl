immutable WeightPointValue{W,P,V}
  weight::W
  point::P
  value::V
end

quadpoints(chart, rule) = quadpoints(chart, rule, Val{dimension(chart)})

function quadpoints(chart, rule, dim::Type{Val{1}})
    u, w = legendre(rule, 0.0, 1.0)
    [meshpoint(chart, u[:,i]) for i in eachindex(w)], w
end

function quadpoints(chart, rule, dim::Type{Val{2}})
    u, w = trgauss(rule)
    [meshpoint(chart, u[:,i]) for i in eachindex(w)], w
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
                qd[i,j][k] = WeightPointValue(w[k],p[k],f(p[k]))
            end
        end
    end

  qd
end
