export patch

function __init__()
    @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin
        @eval import DelimitedFiles
        @eval fn = joinpath(dirname(@__FILE__),"..","examples","cm.csv")
        @eval cm = mapslices(r->(r[1],r[2],r[3]), DelimitedFiles.readdlm(fn, ','), dims=[2])
        @eval function patch(Γ::AbstractMesh, fcr=nothing, caxis=nothing)

            v = vertexarray(Γ)
            c = cellarray(Γ)

            x = v[:,1];    y = v[:,2];    z = v[:,3]
            i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1

            if fcr == nothing
                a = [cartesian(center(chart(Γ,cell)))[3] for cell in cells(Γ)]
            else
                a = fcr
            end

            m, M = extrema(a)
            if caxis != nothing
                m, M = caxis
            end

            if isapprox(m, M)
                n = ones(Integer, size(a))
            else
                n = floor.(Integer, (a.-m)/(M-m)*(length(cm)-1)).+1
                n = clamp.(n, 1, length(cm))
            end
            fc = [cm[i] for i in n]

            s = PlotlyJS.mesh3d(;
                x=x, y=y, z=z,
                i=i, j=j, k=k,
                facecolor=fc,
            )
        end
    end
end
