using CompScienceMeshes, PlotlyJS

Γ = CompScienceMeshes.meshstarpyramid(1.0, 0.5, 5, 0.25, 0.05)

pt1 = PlotlyJS.plot(
    [patch(Γ, opacity=1.0, color="#bcbcbc"), CompScienceMeshes.wireframe(Γ), CompScienceMeshes.normals(Γ)],
    Layout(
        height=400, width=600,
        scene_aspectratio=attr(x=2.0, y=2.0, z=1.0),
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
                eye=attr(x=0, y=1, z=2)
            )
        )
    )
)