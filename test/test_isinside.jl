using Base.Test
using FixedSizeArrays

using CompScienceMeshes

const e0 = Point(0.0, 0.0, 0.0)
const e1 = Point(1.0, 0.0, 0.0)
const e2 = Point(0.0, 1.0, 0.0)
const e3 = Point(0.0, 0.0, 1.0)

p = patch([e0,e1,e2], Val{2})

@test isinside(p, (e0+e1+e2)/3) == true
@test isinside(p, Point(-0.1, 0.3, 0.0)) == false
@test isinside(p, Point(0.3, -0.3, 0.0)) == false
@test isinside(p, Point(0.6, 0.6, 0.0)) == false
@test isinside(p, Point(0.4, 0.6, 0.0)) == false


@test isinclosure(p, (e0+e1+e2)/3) == true
@test isinclosure(p, Point(-0.1, 0.3, 0.0)) == false
@test isinclosure(p, Point(0.3, -0.3, 0.0)) == false
@test isinclosure(p, Point(0.6, 0.6, 0.0)) == false
@test isinclosure(p, Point(0.4, 0.6, 0.0)) == true
