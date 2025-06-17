@testitem "PlotlyJS extension mesh3d/scatter3d" begin
    # using CompScienceMeshes
    using PlotlyBase

    m = meshsphere(radius=1.0, h=0.25)
    t1 = mesh3d(m; opacity=0.5)
    t2 = scatter3d(m; line_width=2, line_color=:black)
    # plt = Plot([t1,t2])
    # display(plt)

    @test t1[:type] == "mesh3d"
    @test t2[:type] == "scatter3d"

end
