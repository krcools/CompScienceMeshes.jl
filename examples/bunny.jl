using CompScienceMeshes
using PlotlyJS

Γ = CompScienceMeshes.meshbunny(0.025; order=3)

# Rotate a quarter turn about the x-axis (ears from +y to +z) and lift onto the xy-plane
# Γ.vertices .= [point(v[1], -v[3], v[2]) for v in Γ.vertices]
# zmin = minimum(v[3] for v in Γ.vertices)
# Γ.vertices .= [v - point(0, 0, zmin) for v in Γ.vertices]

@which CompScienceMeshes.skeleton_fast(Γ, Val{1})
@which CompScienceMeshes.skeleton(Γ, 1)

sk = CompScienceMeshes.skeleton_fast(Γ,Val{1})
simplices = cells(sk)
ctrs = [cartesian(center(chart(sk, cellptr))) for cellptr in sk]
CompScienceMeshes.sort_sfc(ctrs)

@show length(Γ)
Plot(patch(Γ)) |> display
@info "Buidling refinement...."
r = CompScienceMeshes.bisecting_refinement(Γ; order=2)
@info "Refinement built."
@show length(r)
Plot(patch(r)) |> display