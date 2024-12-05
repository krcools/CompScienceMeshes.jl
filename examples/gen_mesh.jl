using CompScienceMeshes, PlotlyJS

Γ = CompScienceMeshes.meshsphere(1.0, 0.3)
Γ = CompScienceMeshes.meshtorus(0.75, 0.25, 0.15)
Γ = CompScienceMeshes.meshstarpyramid(1.0, 0.5, 5, 0.25, 0.05)
Γ = CompScienceMeshes.meshopenbook(π/2, 4, 0.1)

pt1 = PlotlyJS.plot(
    [patch(Γ, opacity=1.0, color="#bcbcbc", lighting = attr(ambient=1)), CompScienceMeshes.wireframe(Γ)],#, CompScienceMeshes.normals(Γ)],
    Layout(
        height=400, width=600,
        scene_aspectratio=attr(x=1.0, y=1.0, z=1.0),
        scene=attr(
            xaxis=attr(
                visible=false,
                showbackground=false
            ),
            yaxis=attr(
                visible=false,
                showbackground=false
            ),
            zaxis=attr(
                visible=false,
                showbackground=false
            ),
            camera=attr(
                up=attr(x=0, y=0, z=1),
                center=attr(x=0, y=0, z=0),
                eye=attr(x=0, y=0.9, z=0.9)
            )
        )
    )
)

PlotlyJS.savefig(pt1, "sphere.pdf", width=600, height=600)