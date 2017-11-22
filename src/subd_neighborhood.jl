export neighborhood,parametric,cartesian,jacobian,normal,shapefuns,shapeders
type subd_point
    cartcoords::Vector{Float64}
    paracoords::Vector{Float64}
    detJ::Float64
    normal::SVector{3,Float64}
    shapes::Matrix{Float64}
    shapeders::Matrix{Float64}
end
Base.getindex(p::subd_point, i::Int) = p.cartcoords[i]
function neighborhood(chart::subd_chart, u)
    vertices = chart.vertices
    Nv = chart.N
    val = Nv - 6
    shapefunc = shape_function(u[1],u[2],val)
    shape_der = shape_function_der(u[1],u[2],val)
    coords = zeros(3)
    t1 = zeros(3)
    t2 = zeros(3)
    for id = 1:3
        for iv = 1:Nv
            coords[id] += vertices[iv][id]*shapefunc[iv]
            t1[id] += vertices[iv][id]*shape_der[iv,1]
            t2[id] += vertices[iv][id]*shape_der[iv,2]
        end
    end
    norm = [t1[2] * t2[3] - t1[3] * t2[2];
            t1[3] * t2[1] - t1[1] * t2[3];
            t1[1] * t2[2] - t1[2] * t2[1]]
    detJ = sqrt(norm[1]*norm[1] + norm[2]*norm[2] + norm[3]*norm[3])
    norm = norm./detJ
    return subd_point(coords,u,detJ,norm,shapefunc,shape_der)
end

function parametric(d::subd_point)
    return d.paracoords
end

function cartesian(d::subd_point)
    v = d.cartcoords
    return v
end

function jacobian(d::subd_point)
    return d.detJ
end

function normal(d::subd_point)
    return d.normal
end

function shapefuns(d::subd_point)
    return d.shapes
end

function shapeder(d::subd_point)
    return d.shapeders
end
