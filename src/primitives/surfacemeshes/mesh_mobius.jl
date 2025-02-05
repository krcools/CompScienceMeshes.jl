"""
    meshmobius(;h)

"""
function meshmobius(;h)

    m = meshrectangle(2pi, 2.0, h)
    m = translate(m, point(-pi, -1, 0))
    for (i,v) in enumerate(m.vertices)
        s,t = v[1], v[2]
        x = 2*point(cos(s), sin(s), 0) + t*(cos(s/2)*point(cos(s), sin(s),0) + sin(s/2)*point(0,0,1))
        m.vertices[i] = x
    end

    n = length(m)
    m2 = deepcopy(m)
    m3 = weld(m, m2; glueop=identity)
    return Mesh(m3.vertices, m3.faces[n+1:end])
end