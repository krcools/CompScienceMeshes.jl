module PlotlyBaseExt

using CompScienceMeshes
using PlotlyBase
# using PlotlyJS.PlotlyBase

function CompScienceMeshes.patch(Γ::CompScienceMeshes.AbstractMesh, fcr=nothing;
    caxis=nothing, showscale=true, color="red", kwargs...)

    v = vertexarray(Γ)
    c = cellarray(Γ)

    x = v[:,1];    y = v[:,2];    z = v[:,3]
    i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1


    if fcr != nothing
        m, M = extrema(fcr)
        if caxis != nothing
            m, M = caxis
        end

        s = PlotlyBase.mesh3d(;
            x=x, y=y, z=z,
            i=i, j=j, k=k,    
            intensitymode="cell",
            intensity=fcr,
            colorscale="Viridis",
            showscale=showscale,
            cmin=m,
            cmax=M,
            kwargs...
        )
    else
        s = PlotlyBase.mesh3d(;
            x=x, y=y, z=z,
            i=i, j=j, k=k,
            color=color,
            kwargs...
        )
    end
    return s
end

function CompScienceMeshes.patch(a::Vector{<:CompScienceMeshes.Simplex}; kwargs...)
    vertices = reduce(vcat, [v.vertices for v in a])
    faces = collect(SVector(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3) for i in 1:length(a))
    mesh = CompScienceMeshes.Mesh(vertices, faces)
    return CompScienceMeshes.patch(mesh; kwargs...)
end

function PlotlyBase.mesh3d(geo::CompScienceMeshes.AbstractMesh; kwargs...)

    v = CompScienceMeshes.vertexarray(geo)
    c = CompScienceMeshes.cellarray(geo)

    x = v[:,1];    y = v[:,2];    z = v[:,3]
    i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1

    return PlotlyBase.mesh3d(;
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        kwargs...)
end

function PlotlyBase.cone(mesh, arrows; sizeref=2, kwargs...)
    centers = [cartesian(center(chart(mesh,cell))) for cell in mesh]
    x = getindex.(centers,1)
    y = getindex.(centers,2)
    z = getindex.(centers,3)
    u = getindex.(arrows,1)
    v = getindex.(arrows,2)
    w = getindex.(arrows,3)
    PlotlyBase.cone(x=x,y=y,z=z,u=u,v=v,w=w,sizemode="absolute", sizeref=sizeref; kwargs...)
end

function CompScienceMeshes.normals(mesh; kwargs...)
    nrmls = [normal(chart(mesh,cell)) for cell in mesh]
    return PlotlyBase.cones(mesh, nrmls; kwargs...)
end

function CompScienceMeshes.wireframe(edges; width=1, color="rgb(0,0,0)")
    edges = skeleton(edges,1)
    T = coordtype(edges)
    x = T[]
    y = T[]
    z = T[]
    for edge in edges
        chrt = chart(edges, edge)
        v1, v2 = chrt.vertices
        append!(x, [v1[1],v2[1],NaN])
        append!(y, [v1[2],v2[2],NaN])
        append!(z, [v1[3],v2[3],NaN])
    end
    return PlotlyBase.scatter3d(x=x,y=y,z=z,mode="lines",
        line=PlotlyBase.attr(
            color=color,
            width=width
        )
    )
end

function PlotlyBase.scatter3d(edges::CompScienceMeshes.AbstractMesh{3,2}; kwargs...)
    T = CompScienceMeshes.coordtype(edges)
    x = T[]
    y = T[]
    z = T[]
    for edge in edges
        chrt = CompScienceMeshes.chart(edges, edge)
        v1, v2 = chrt.vertices
        append!(x, [v1[1],v2[1],NaN])
        append!(y, [v1[2],v2[2],NaN])
        append!(z, [v1[3],v2[3],NaN])
    end
    return PlotlyBase.scatter3d(x=x, y=y, z=z, mode="lines"; kwargs...)
end

function PlotlyBase.scatter3d(surf::CompScienceMeshes.AbstractMesh{3,3}; kwargs...)
    PlotlyBase.scatter3d(CompScienceMeshes.skeleton(surf, 1); kwargs...)
end

end