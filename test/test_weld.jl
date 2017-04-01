using Base.Test

using StaticArrays
using CompScienceMeshes

m1 = meshrectangle(1.0, 1.0, 0.5)
m2 = meshrectangle(1.0, 1.0, 0.5)
translate!(m2, point(0,1,0))
#m2.vertices += point(0.0, 1.0, 0.0)

m = weld(m1, m2)

@test numcells(skeleton(m,0)) == 15


## The mysterious case of the cuboid
using CompScienceMeshes
using Base.Test

width, height, depth = 1.0, 1.0, 1.0

h = 1/4

m = meshrectangle(width, height, h, 3)

m_front = rotate(m, 0.5*pi*[1,0,0])
m_left  = rotate(m, -0.5*pi*[0,1,0])

m_back  = translate(m_front, point(0,height,0))
fliporientation!(m_back)
m_right = translate(m_left, point(width,0,0))
fliporientation!(m_right)

m_bottom = fliporientation(m)

m_fine = meshrectangle(width, height, h)
m_top = translate(m_fine, point(0,0,height))
m_open = weld(m_front, m_left, m_back, m_right, m_bottom)

#include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))
#patch(m_top)
ot = weld(m_open, m_top)
#patch(ot)

nvo = numcells(skeleton(m_open,0))
nvt = numcells(skeleton(m_top,0))
nvm = numcells(skeleton(ot,0))
nvb = numcells(skeleton(boundary(m_top),0))
@test nvm == nvo + nvt - nvb
