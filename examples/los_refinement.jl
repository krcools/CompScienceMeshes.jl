using CompScienceMeshes

coarse = meshsphere(1.0, 0.45)
fine = CompScienceMeshes.lineofsight_refinement(coarse)
vtoc, nc = CompScienceMeshes.vertextocellmap(fine)

colors = zeros(numcells(fine))
for i in 1:numvertices(coarse)
    c = rand()
    for j in 1 : (nc[i])
        colors[vtoc[i,j]] = c
    end
end

import Plotly
Plotly.plot(patch(fine.mesh, colors))
