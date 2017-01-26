using CompScienceMeshes

Γ = meshsphere(1.0, 0.3)
γ = bisecting_refinement(Γ)

include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
p = patch(γ)
PlotlyJS.plot(p)
