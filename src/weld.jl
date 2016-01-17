export weld

using CollisionDetection
using CompScienceMeshes

function weld(Γ₁, Γ₂)

  tree = CollisionDetection.octree(Γ₁.vertices)

  idmap = collect(numvertices(Γ₁) + (1:numvertices(Γ₂)))
  for (j,v) in enumerate(Γ₂.vertices)
    I = find(tree, v)
    isempty(I) || (idmap[j] = I[1])
  end

  V = [Γ₁.vertices; Γ₂.vertices]
  F = [Γ₁.faces; Γ₂.faces]
  for i in 1:numcells(Γ₂)
    c = Γ₂.faces[i]
    F[numcells(Γ₁) + i] = Vec(idmap[c[1]], idmap[c[2]], idmap[c[3]])
  end

  return Mesh(V,F)

end
