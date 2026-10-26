using CompScienceMeshes
using PlotlyBase

m = meshsphere(radius=1.0, h=0.3)
n = collect(normal(chart(m,c)) for c in m)

tr1 = mesh3d(m; opacity=0.5)
tr2 = cone(m, n; sizeref=2)
tr3 = scatter3d(skeleton(m,1); line_width=4, line_color="black")
Plot([tr1,tr2,tr3])
