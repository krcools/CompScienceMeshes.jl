import PlotlyJS
import Colors

const cm = mapslices(
    r->Colors.RGB(r[1],r[2],r[3]),
    readcsv(Pkg.dir("CompScienceMeshes","examples","cm.csv")),
    [2])

function patch(Γ, fcr)

    v = vertexarray(Γ)
    c = cellarray(Γ)

    x = v[:,1]; y = v[:,2]; z = v[:,3]
    i = c[:,1]-1; j = c[:,2]-1; k = c[:,3]-1

    a = Float64[sqrt(real(dot(f,f))) for f in fcr]
    m, M = extrema(a)

    n = floor(Integer, (a-m)/(M-m)*(length(cm)-1))+1
    fc = [cm[i] for i in n]

    s = PlotlyJS.mesh3d(;
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        facecolor=fc,
    )

    #PlotlyJS.plot([s])
    s
end
