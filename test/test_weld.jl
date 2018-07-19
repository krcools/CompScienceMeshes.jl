using Test

using StaticArrays
using CompScienceMeshes

#include(Pkg.dir("CompScienceMeshes","examples","matlab_patches.jl"))

m1 = meshrectangle(1.0, 1.0, 0.5)
m2 = meshrectangle(1.0, 1.0, 0.5)
translate!(m2, point(0,1,0))
#m2.vertices += point(0.0, 1.0, 0.0)

m = weld(m1, m2)

@test numcells(skeleton(m,0)) == 15


## The mysterious case of the cuboid
using CompScienceMeshes
using Test

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

#patch(m_top)
ot = weld(m_open, m_top)
#patch(ot)

nvo = numcells(skeleton(m_open,0))
nvt = numcells(skeleton(m_top,0))
nvm = numcells(skeleton(ot,0))
nvb = numcells(skeleton(boundary(m_top),0))
@test nvm == nvo + nvt - nvb


# test that consecutive welding happens correctly
function meshopenbox(a, h)

    base = meshrectangle(a, a, h, 3)
    front = rotate(base, point(π/2,0,0))
    back = -translate(front, point(0,a,0))
    left = rotate(base, point(0,-π/2,0))
    right = -translate(left, point(a,0,0))

    return weld(-base, front, right, back, left)
end

function meshuboidwithdivision(a, h)

    m1 = meshopenbox(a,h)
    translate!(m1, point(0,0,-a))

    m2 = deepcopy(m1)
    m2.vertices = [point(v[1],v[2],-v[3]) for v in m2.vertices]

    m3 = meshrectangle(a, a, h, 3)

    m = weld(m1,m2,m3)
    return m
end

a = 1.0
m = meshopenbox(a,a)
@test numcells(skeleton(m,0)) == 8

m = meshuboidwithdivision(a,a)
verts = skeleton(m, 0)
@test numcells(verts) == 12

##
a = h = 1.0
base = meshrectangle(a, a, h, 3)
front = rotate(base, point(π/2,0,0))
back = -translate(front, point(0,a,0))
left = rotate(base, point(0,-π/2,0))
right = -translate(left, point(a,0,0))

m1= weld(-base, front)
@test numcells(skeleton(m1,0)) == 6
m2 = weld(m1, right)
@test numcells(skeleton(m2,0)) == 7

m3 = weld(front, right)
@test numcells(skeleton(m3,0)) == 6


#patch(m)
