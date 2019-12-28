using CompScienceMeshes

m = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
@show numcells(m)

br = barycentric_refinement(m)
@show numcells(br)

@assert numcells(br) == 24*numcells(m)

V = 155
sm = submesh(tet -> V in tet, br.mesh)
@show numcells(sm)

bndry = boundary(sm)
import PlotlyJS
PlotlyJS.plot(patch(bndry))

E = 100
Edges = skeleton(m,1)
sm = submesh(br.mesh) do tet
    cells(Edges)[E][1] in tet && return true
    cells(Edges)[E][2] in tet && return true
    return false
end

PlotlyJS.plot(patch(boundary(sm)))
