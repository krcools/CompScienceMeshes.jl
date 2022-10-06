using CompScienceMeshes

m = CompScienceMeshes.tetmeshsphere(1.0, 0.35)
@show length(m)

br = barycentric_refinement(m)
@show length(br)

@assert length(br) == 24*length(m)

V = 155
# sm = submesh(tet -> V in tet, br.mesh)
sm = submesh(br.mesh) do m,p
    inds = CompScienceMeshes.indices(m,p)
    return V in inds
end
@show length(sm)

bndry = boundary(sm)
import Plotly
Plotly.plot(patch(bndry))

E = 100
Edges = skeleton(m,1)
sm = submesh(br.mesh) do m,p
    tet = CompScienceMeshes.indices(m,p)
    edge = CompScienceMeshes.indices(Edges,E)
    edge[1] in tet && return true
    edge[2] in tet && return true
    return false
end

bnd = boundary(sm)
Plotly.plot(patch(bnd, range(0,1,length=length(bnd))))
