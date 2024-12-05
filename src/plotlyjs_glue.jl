export patch
export normalcones
export wireframe

function __init__()
    @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin

        @eval function patch(Γ::CompScienceMeshes.AbstractMesh, fcr=nothing;
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
        
                s = PlotlyJS.mesh3d(;
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
                s = PlotlyJS.mesh3d(;
                    x=x, y=y, z=z,
                    i=i, j=j, k=k,
                    color=color,
                    kwargs...
                )
            end
            return s
        end

        @eval function patch(a::Vector{<:Simplex}; kwargs...)
            vertices = reduce(vcat, [v.vertices for v in a])
            faces = collect(SVector(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3) for i in 1:length(a))
            # faces = reduce(vcat, SVector(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3) for i in 1:length(a))
            mesh = Mesh(vertices, faces)
            return patch(mesh; kwargs...)
        end

        @eval function cones(mesh, arrows; sizeref=2, kwargs...)
            # normals = [normal(chart(mesh,cell)) for cell in mesh]
            centers = [cartesian(center(chart(mesh,cell))) for cell in mesh]
            x = getindex.(centers,1)
            y = getindex.(centers,2)
            z = getindex.(centers,3)
            u = getindex.(arrows,1)
            v = getindex.(arrows,2)
            w = getindex.(arrows,3)
            PlotlyJS.cone(x=x,y=y,z=z,u=u,v=v,w=w,sizemode="absolute", sizeref=sizeref; kwargs...)
        end

        @eval function normals(mesh; kwargs...)
            nrmls = [normal(chart(mesh,cell)) for cell in mesh]
            return CompScienceMeshes.cones(mesh, nrmls; kwargs...)
        end

        @eval function wireframe(edges; width=1, color="rgb(0,0,0)")
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
            return PlotlyJS.scatter3d(x=x,y=y,z=z,mode="lines",
                line=PlotlyJS.attr(
                    color=color,
                    width=width
                )
            )
        end
    end
end