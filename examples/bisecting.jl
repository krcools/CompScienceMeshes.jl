using CompScienceMeshes

Γ = meshsphere(radius=1.0, h=0.3)
γ = bisecting_refinement(Γ)

import Pkg
# include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
import Plotly
p = patch(γ; opacity=0.9)
w0 = CompScienceMeshes.wireframe(skeleton(Γ,1); width=4.0)
w1 = CompScienceMeshes.wireframe(skeleton(γ,1); width=1.0)
Plotly.plot([p,w0,w1])
