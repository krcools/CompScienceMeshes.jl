using Base.Test
using FixedSizeArrays

import CompScienceMeshes
CM = CompScienceMeshes

m1 = CM.meshrectangle(1.0, 1.0, 0.5)
m2 = CM.meshrectangle(1.0, 1.0, 0.5)
m2.vertices += Point(0.0, 1.0, 0.0)

m = CM.weld(m1, m2)

@test length(CM.cells(m,0)) == 15
