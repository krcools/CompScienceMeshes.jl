using Base.Test

using CompScienceMeshes

m1 = meshrectangle(1.0, 1.0, 0.5)
m2 = meshrectangle(1.0, 1.0, 0.5)
translate!(m2, point(0,1,0))
#m2.vertices += point(0.0, 1.0, 0.0)

m = weld(m1, m2)

@test numcells(skeleton(m,0)) == 15
