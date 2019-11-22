using CompScienceMeshes
import PlotlyJS

Faces = meshsphere(1.0, 0.55)

faces = barycentric_refinement(Faces)
PlotlyJS.plot(patch(faces.mesh))
# error("top")
NV = numvertices(Faces)

# Nodes = skeleton(Faces, 0)
# node_ctrs = [vertices(Nodes)[node][1] for node in cells(Nodes)]
# # C = [c[i] for c in node_ctrs, i in 1:3]
# Nodes.faces = Nodes.faces[CompScienceMeshes.sort_sfc(node_ctrs)]
#
# fine = Mesh(verts, fcs)
# D = connectivity(Nodes, fine)
# rows, vals = rowvals(D), nonzeros(D)
# sorted_fcs = Vector{celltype(mesh)}()
# for (i,Node) in enumerate(cells(Nodes))
#     for k in nzrange(D,i)
#         j = rows[k]
#         push!(sorted_fcs, fcs[j])
#     end
# end



vertices(Faces)[1:NV] == vertices(faces)[1:NV]

Nodes = skeleton(Faces,0)

D = connectivity(Nodes, faces)

import PlotlyJS
numcells(faces)
nf = 224
t1 = patch(Mesh(faces.mesh.vertices, faces.mesh.faces[1:nf]), rand(numcells(faces)))
PlotlyJS.plot(t1)

ctrs = [cartesian(center(chart(Nodes,c))) for c in cells(Nodes)]
Ctrs = [ctr[i] for ctr in ctrs, i in 1:3]
nc = div(length(ctrs),3)
t2 = PlotlyJS.scatter3d(x=Ctrs[1:nc,1], y=Ctrs[1:nc,2], z=Ctrs[1:nc,3])

Edges = skeleton(Faces,1)
ctrs = [cartesian(center(chart(Edges,c))) for c in cells(Edges)]
Ctrs = [ctr[i] for ctr in ctrs, i in 1:3]
nc = div(length(ctrs),3)
t3 = PlotlyJS.scatter3d(x=Ctrs[1:nc,1], y=Ctrs[1:nc,2], z=Ctrs[1:nc,3])

PlotlyJS.plot([t1,t2,t3])

aux = skeleton(faces,2)

t4 = patch(Mesh(aux.vertices, aux.faces[1:nf]), rand(nf))

PlotlyJS.plot(t4)
