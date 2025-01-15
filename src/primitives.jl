"""
    meshrectangle(width, height, delta, udim)

Create a mesh for a rectangle of width (along the x-axis) `width` and height (along
    the y-axis) `height`.

The target edge size is `delta` and the dimension of the
    embedding universe is `udim` (>= 2).

The mesh is oriented such that the normal is pointing down. This is subject to change.
"""
function meshrectangle(width::T, height::T, delta::T, udim=3; structured=true) where T
  if structured == false
	  @assert udim==3 "Only 3D Unstructured mesh currently supported"
	  return meshrectangle_unstructured(width, height, delta)
  end

  if structured == :quadrilateral || structured == "quadrilateral"
    @assert udim==3 "Only 3D Unstructured mesh currently supported"
    return meshrectangle_quadrilateral(width, height, delta)
  end

  PT = SVector{udim,T}
  CT = SVector{3,Int}

    @assert 2 <= udim

    nx = round(Int, ceil(width/delta));  nx = max(nx,1); dx = width/nx
    ny = round(Int, ceil(height/delta)); ny = max(ny,1); dy = height/ny

    xs = (0:nx) * dx
    ys = (0:ny) * dy

    vertices = zeros(PT, (nx+1)*(ny+1))
    k = 1
    for x in xs
        for y in ys
            p = zeros(T, udim)
            p[1] = x
            p[2] = y
            vertices[k] = PT(p)
            k += 1
        end
    end

    faces = zeros(CT, 2*nx*ny)
    k = 1
    for i in 1 : nx
        for j in 1 : ny
            v11 = (i-1)*(ny+1) + j
            v12 = i*(ny+1) + j
            v21 = (i-1)*(ny+1) + (j+1)
            v22 = i*(ny+1) + (j+1)
            faces[k]   = CT(v11,v21,v12)
            faces[k+1] = CT(v21,v22,v12)
            k += 2
        end
    end

    Mesh(vertices, faces)
end

function meshrectangle_quadrilateral(width::T, height::T, delta::T) where {T}

    nx = round(Int, ceil(width/delta));  nx = max(nx,1); dx = width/nx
    ny = round(Int, ceil(height/delta)); ny = max(ny,1); dy = height/ny

    xs = range(0, width, step=dx)
    ys = range(0, height, step=dy)

    vertices = Vector{SVector{3,T}}(undef, (nx+1)*(ny+1))
    k = 1
    for i in 1:nx+1
        for j in 1:ny+1
            vertices[k] = SVector{3,T}((i-1)*dx, (j-1)*dy, 0)
            k += 1
    end end

    faces = Vector{SVector{4,Int}}(undef, nx*ny)
    k = 1
    for i in 1:nx
        for j in 1:ny
            faces[k] = SVector{4,Int}(
                (i-1)*(ny+1) + j, i*(ny+1) + j, i*(ny+1) + j+1, (i-1)*(ny+1) + j+1)
            k += 1
    end end

    return QuadMesh(vertices, faces)
end

@testitem "QuadMesh rectangle" begin
    m = meshrectangle(2.0, 2.0, 1.0; structured=:quadrilateral)
    @test length(m) == 4
    ch = chart(m, first(m))
    p = neighborhood(ch, point(0.5, 0.5))
    @test cartesian(p) ≈ point(0.5, 0.5, 0.0)
end

"""
    meshtorus(majorradius, minorradius, h)

Create a mesh of a torus of 2 radii `majorradius` and `minorradius`. The axis of symmetry is ̂z.

The target edge size is `h`.
"""
function meshtorus(majorradius, minorradius, h)
    @assert minorradius < majorradius
    
    fno = tempname() * ".msh"

    gmsh.initialize()
    gmsh.model.add("torus")
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 1)
    gmsh.model.geo.addPoint(0.0, 0.0, minorradius, h, 2)
    gmsh.model.geo.addPoint(0.0, 0.0, -minorradius, h, 3)

    gmsh.model.geo.addPoint(majorradius, 0.0, 0.0, h, 4)
    gmsh.model.geo.addPoint(majorradius + minorradius, 0.0, 0.0, h, 5)
    gmsh.model.geo.addPoint(majorradius, 0.0, minorradius, h, 6)
    gmsh.model.geo.addPoint(majorradius - minorradius, 0.0, 0.0, h, 7)
    gmsh.model.geo.addPoint(majorradius, 0.0, -minorradius, h, 8)
    gmsh.model.geo.addCircleArc(5, 4, 6, 1)
    gmsh.model.geo.addCircleArc(6, 4, 7, 2)
    gmsh.model.geo.addCircleArc(7, 4, 8, 3)
    gmsh.model.geo.addCircleArc(8, 4, 5, 4)

    gmsh.model.geo.addPoint(0.0, majorradius, 0.0, h, 9)
    gmsh.model.geo.addPoint(0.0, majorradius + minorradius, 0.0, h, 10)
    gmsh.model.geo.addPoint(0.0, majorradius, minorradius, h, 11)
    gmsh.model.geo.addPoint(0.0, majorradius - minorradius, 0.0, h, 12)
    gmsh.model.geo.addPoint(0.0, majorradius, -minorradius, h, 13)
    gmsh.model.geo.addCircleArc(10, 9, 11, 5)
    gmsh.model.geo.addCircleArc(11, 9, 12, 6)
    gmsh.model.geo.addCircleArc(12, 9, 13, 7)
    gmsh.model.geo.addCircleArc(13, 9, 10, 8)

    gmsh.model.geo.addPoint(-majorradius, 0.0, 0.0, h, 14)
    gmsh.model.geo.addPoint(-majorradius - minorradius, 0.0, 0.0, h, 15)
    gmsh.model.geo.addPoint(-majorradius, 0.0, minorradius, h, 16)
    gmsh.model.geo.addPoint(-majorradius + minorradius, 0.0, 0.0, h, 17)
    gmsh.model.geo.addPoint(-majorradius, 0.0, -minorradius, h, 18)
    gmsh.model.geo.addCircleArc(15, 14, 16, 9)
    gmsh.model.geo.addCircleArc(16, 14, 17, 10)
    gmsh.model.geo.addCircleArc(17, 14, 18, 11)
    gmsh.model.geo.addCircleArc(18, 14, 15, 12)

    gmsh.model.geo.addPoint(0.0, -majorradius, 0.0, h, 19)
    gmsh.model.geo.addPoint(0.0, -majorradius - minorradius, 0.0, h, 20)
    gmsh.model.geo.addPoint(0.0, -majorradius, minorradius, h, 21)
    gmsh.model.geo.addPoint(0.0, -majorradius + minorradius, 0.0, h, 22)
    gmsh.model.geo.addPoint(0.0, -majorradius, -minorradius, h, 23)
    gmsh.model.geo.addCircleArc(20, 19, 21, 13)
    gmsh.model.geo.addCircleArc(21, 19, 22, 14)
    gmsh.model.geo.addCircleArc(22, 19, 23, 15)
    gmsh.model.geo.addCircleArc(23, 19, 20, 16)
    
    gmsh.model.geo.addCircleArc(5, 1, 10, 17)
    gmsh.model.geo.addCircleArc(6, 2, 11, 18)
    gmsh.model.geo.addCircleArc(7, 1, 12, 19)
    gmsh.model.geo.addCircleArc(8, 3, 13, 20)

    gmsh.model.geo.addCircleArc(10, 1, 15, 21)
    gmsh.model.geo.addCircleArc(11, 2, 16, 22)
    gmsh.model.geo.addCircleArc(12, 1, 17, 23)
    gmsh.model.geo.addCircleArc(13, 3, 18, 24)

    gmsh.model.geo.addCircleArc(15, 1, 20, 25)
    gmsh.model.geo.addCircleArc(16, 2, 21, 26)
    gmsh.model.geo.addCircleArc(17, 1, 22, 27)
    gmsh.model.geo.addCircleArc(18, 3, 23, 28)

    gmsh.model.geo.addCircleArc(20, 1, 5, 29)
    gmsh.model.geo.addCircleArc(21, 2, 6, 30)
    gmsh.model.geo.addCircleArc(22, 1, 7, 31)
    gmsh.model.geo.addCircleArc(23, 3, 8, 32)

    gmsh.model.geo.addCurveLoop([-1, -18, 5, 17], 33)
    gmsh.model.geo.addSurfaceFilling([33], 1)
    gmsh.model.geo.addCurveLoop([-2, -19, 6, 18], 34)
    gmsh.model.geo.addSurfaceFilling([34], 2)
    gmsh.model.geo.addCurveLoop([-3, -20, 7, 19], 35)
    gmsh.model.geo.addSurfaceFilling([35], 3)
    gmsh.model.geo.addCurveLoop([-4, -17, 8, 20], 36)
    gmsh.model.geo.addSurfaceFilling([36], 4)

    gmsh.model.geo.addCurveLoop([-5, -22, 9, 21], 37)
    gmsh.model.geo.addSurfaceFilling([37], 5)
    gmsh.model.geo.addCurveLoop([-6, -23, 10, 22], 38)
    gmsh.model.geo.addSurfaceFilling([38], 6)
    gmsh.model.geo.addCurveLoop([-7, -24, 11, 23], 39)
    gmsh.model.geo.addSurfaceFilling([39], 7)
    gmsh.model.geo.addCurveLoop([-8, -21, 12, 24], 40)
    gmsh.model.geo.addSurfaceFilling([40], 8)
    
    gmsh.model.geo.addCurveLoop([-9, -26, 13, 25], 41)
    gmsh.model.geo.addSurfaceFilling([41], 9)
    gmsh.model.geo.addCurveLoop([-10, -27, 14, 26], 42)
    gmsh.model.geo.addSurfaceFilling([42], 10)
    gmsh.model.geo.addCurveLoop([-11, -28, 15, 27], 43)
    gmsh.model.geo.addSurfaceFilling([43], 11)
    gmsh.model.geo.addCurveLoop([-12, -25, 16, 28], 44)
    gmsh.model.geo.addSurfaceFilling([44], 12)

    gmsh.model.geo.addCurveLoop([-13, -30, 1, 29], 45)
    gmsh.model.geo.addSurfaceFilling([45], 13)
    gmsh.model.geo.addCurveLoop([-14, -31, 2, 30], 46)
    gmsh.model.geo.addSurfaceFilling([46], 14)
    gmsh.model.geo.addCurveLoop([-15, -32, 3, 31], 47)
    gmsh.model.geo.addSurfaceFilling([47], 15)
    gmsh.model.geo.addCurveLoop([-16, -29, 4, 32], 48)
    gmsh.model.geo.addSurfaceFilling([48], 16)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end




"""
    meshcuboid1hole(width, height, holewidth, h)

Create a mesh of a cuboid of size `width × width × height` with a hole of size `holewidth × holewidth × height` (̂x × ̂y × ̂z, resp.).

The target edge size is `h`.
"""
function meshcuboid1hole(width, height, holewidth, h)
    @assert holewidth < width
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("squaretorus")

    # bottom plate
    gmsh.model.geo.addPoint(holewidth/2, -holewidth/2, -height/2, h, 1)
    gmsh.model.geo.addPoint(holewidth/2, holewidth/2, -height/2, h, 2)
    gmsh.model.geo.addPoint(-holewidth/2, holewidth/2, -height/2, h, 3)
    gmsh.model.geo.addPoint(-holewidth/2, -holewidth/2, -height/2, h, 4)
    gmsh.model.geo.addPoint(width/2, -width/2, -height/2, h, 5)
    gmsh.model.geo.addPoint(width/2, width/2, -height/2, h, 6)
    gmsh.model.geo.addPoint(-width/2, width/2, -height/2, h, 7)
    gmsh.model.geo.addPoint(-width/2, -width/2, -height/2, h, 8)

    gmsh.model.geo.addLine(2, 3, 1)
    gmsh.model.geo.addLine(3, 4, 2)
    gmsh.model.geo.addLine(4, 1, 3)
    gmsh.model.geo.addLine(1, 2, 4)
    gmsh.model.geo.addCurveLoop([-1, -2, -3, -4], 101)

    gmsh.model.geo.addLine(6, 7, 5)
    gmsh.model.geo.addLine(7, 8, 6)
    gmsh.model.geo.addLine(8, 5, 7)
    gmsh.model.geo.addLine(5, 6, 8)
    gmsh.model.geo.addCurveLoop([-5, -6, -7, -8], 102)
    gmsh.model.geo.addPlaneSurface([-101, -102], 1)

    # top plate
    gmsh.model.geo.addPoint(holewidth/2, -holewidth/2, height/2, h, 9)
    gmsh.model.geo.addPoint(holewidth/2, holewidth/2, height/2, h, 10)
    gmsh.model.geo.addPoint(-holewidth/2, holewidth/2, height/2, h, 11)
    gmsh.model.geo.addPoint(-holewidth/2, -holewidth/2, height/2, h, 12)
    gmsh.model.geo.addPoint(width/2, -width/2, height/2, h, 13)
    gmsh.model.geo.addPoint(width/2, width/2, height/2, h, 14)
    gmsh.model.geo.addPoint(-width/2, width/2, height/2, h, 15)
    gmsh.model.geo.addPoint(-width/2, -width/2, height/2, h, 16)

    gmsh.model.geo.addLine(10, 11, 9)
    gmsh.model.geo.addLine(11, 12, 10)
    gmsh.model.geo.addLine(12, 9, 11)
    gmsh.model.geo.addLine(9, 10, 12)
    gmsh.model.geo.addCurveLoop([9, 10, 11, 12], 103)
    
    gmsh.model.geo.addLine(14, 15, 13)
    gmsh.model.geo.addLine(15, 16, 14)
    gmsh.model.geo.addLine(16, 13, 15)
    gmsh.model.geo.addLine(13, 14, 16)
    gmsh.model.geo.addCurveLoop([13, 14, 15, 16], 104)
    gmsh.model.geo.addPlaneSurface([-103, -104], 2)

    # sides
    gmsh.model.geo.addLine(2, 10, 17)
    gmsh.model.geo.addLine(3, 11, 18)
    gmsh.model.geo.addLine(4, 12, 19)
    gmsh.model.geo.addLine(1, 9, 20)
    gmsh.model.geo.addLine(6, 14, 21)
    gmsh.model.geo.addLine(7, 15, 22)
    gmsh.model.geo.addLine(8, 16, 23)
    gmsh.model.geo.addLine(5, 13, 24)

    gmsh.model.geo.addCurveLoop([-12, 17, 4, -20], 105)
    gmsh.model.geo.addPlaneSurface([-105], 3)

    gmsh.model.geo.addCurveLoop([1, 18, -9, -17], 106)
    gmsh.model.geo.addPlaneSurface([-106], 4)

    gmsh.model.geo.addCurveLoop([-18, -10, 19, 2], 107)
    gmsh.model.geo.addPlaneSurface([-107], 5)

    gmsh.model.geo.addCurveLoop([-19, -11, 20, 3], 108)
    gmsh.model.geo.addPlaneSurface([-108], 6)

    gmsh.model.geo.addCurveLoop([24, 16, -21, -8], 109)
    gmsh.model.geo.addPlaneSurface([-109], 7)

    gmsh.model.geo.addCurveLoop([-5, -22, 13, 21], 110)
    gmsh.model.geo.addPlaneSurface([-110], 8)

    gmsh.model.geo.addCurveLoop([-6, -23, 14, 22], 111)
    gmsh.model.geo.addPlaneSurface([-111], 9)

    gmsh.model.geo.addCurveLoop([23, 15, -24, -7], 112)
    gmsh.model.geo.addPlaneSurface([-112], 10)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end



"""
    meshcuboid4holes(width, height, holewidth, h)

Create a mesh of a cuboid of size `width × width × height` with 4 holes of size `holewidth × holewidth × height` (̂x × ̂y × ̂z, resp.).

The target edge size is `h`.
"""

function meshcuboid4holes(width, height, holewidth, h)
    @assert 2*holewidth < width
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("squaretorus4holes")

    # bottom plate
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 - holewidth/2, -height/2, h, 1)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 + holewidth/2, -height/2, h, 2)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 + holewidth/2, -height/2, h, 3)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 - holewidth/2, -height/2, h, 4)

    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 - holewidth/2, -height/2, h, 5)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 + holewidth/2, -height/2, h, 6)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 + holewidth/2, -height/2, h, 7)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 - holewidth/2, -height/2, h, 8)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 - holewidth/2, -height/2, h, 9)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 + holewidth/2, -height/2, h, 10)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 + holewidth/2, -height/2, h, 11)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 - holewidth/2, -height/2, h, 12)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 - holewidth/2, -height/2, h, 13)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 + holewidth/2, -height/2, h, 14)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 + holewidth/2, -height/2, h, 15)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 - holewidth/2, -height/2, h, 16)

    gmsh.model.geo.addPoint(width/2, -width/2, -height/2, h, 17)
    gmsh.model.geo.addPoint(width/2, width/2, -height/2, h, 18)
    gmsh.model.geo.addPoint(-width/2, width/2, -height/2, h, 19)
    gmsh.model.geo.addPoint(-width/2, -width/2, -height/2, h, 20)

    gmsh.model.geo.addLine(2, 3, 1)
    gmsh.model.geo.addLine(3, 4, 2)
    gmsh.model.geo.addLine(4, 1, 3)
    gmsh.model.geo.addLine(1, 2, 4)
    gmsh.model.geo.addCurveLoop([-1, -2, -3, -4], 101)

    gmsh.model.geo.addLine(6, 7, 5)
    gmsh.model.geo.addLine(7, 8, 6)
    gmsh.model.geo.addLine(8, 5, 7)
    gmsh.model.geo.addLine(5, 6, 8)
    gmsh.model.geo.addCurveLoop([-5, -6, -7, -8], 102)

    gmsh.model.geo.addLine(10, 11, 9)
    gmsh.model.geo.addLine(11, 12, 10)
    gmsh.model.geo.addLine(12, 9, 11)
    gmsh.model.geo.addLine(9, 10, 12)
    gmsh.model.geo.addCurveLoop([-9, -10, -11, -12], 103)

    gmsh.model.geo.addLine(14, 15, 13)
    gmsh.model.geo.addLine(15, 16, 14)
    gmsh.model.geo.addLine(16, 13, 15)
    gmsh.model.geo.addLine(13, 14, 16)
    gmsh.model.geo.addCurveLoop([-13, -14, -15, -16], 104)

    gmsh.model.geo.addLine(18, 19, 17)
    gmsh.model.geo.addLine(19, 20, 18)
    gmsh.model.geo.addLine(20, 17, 19)
    gmsh.model.geo.addLine(17, 18, 20)
    gmsh.model.geo.addCurveLoop([-17, -18, -19, -20], 105)

    gmsh.model.geo.addPlaneSurface([-101, -102, -103, -104, -105], 1)

    # top plate
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 - holewidth/2, height/2, h, 21)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 + holewidth/2, height/2, h, 22)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 + holewidth/2, height/2, h, 23)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 - holewidth/2, height/2, h, 24)

    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 - holewidth/2, height/2, h, 25)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 + holewidth/2, height/2, h, 26)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 + holewidth/2, height/2, h, 27)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 - holewidth/2, height/2, h, 28)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 - holewidth/2, height/2, h, 29)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 + holewidth/2, height/2, h, 30)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 + holewidth/2, height/2, h, 31)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 - holewidth/2, height/2, h, 32)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 - holewidth/2, height/2, h, 33)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 + holewidth/2, height/2, h, 34)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 + holewidth/2, height/2, h, 35)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 - holewidth/2, height/2, h, 36)

    gmsh.model.geo.addPoint(width/2, -width/2, height/2, h, 37)
    gmsh.model.geo.addPoint(width/2, width/2, height/2, h, 38)
    gmsh.model.geo.addPoint(-width/2, width/2, height/2, h, 39)
    gmsh.model.geo.addPoint(-width/2, -width/2, height/2, h, 40)

    gmsh.model.geo.addLine(22, 23, 21)
    gmsh.model.geo.addLine(23, 24, 22)
    gmsh.model.geo.addLine(24, 21, 23)
    gmsh.model.geo.addLine(21, 22, 24)
    gmsh.model.geo.addCurveLoop([-21, -22, -23, -24], 106)

    gmsh.model.geo.addLine(26, 27, 25)
    gmsh.model.geo.addLine(27, 28, 26)
    gmsh.model.geo.addLine(28, 25, 27)
    gmsh.model.geo.addLine(25, 26, 28)
    gmsh.model.geo.addCurveLoop([-25, -26, -27, -28], 107)

    gmsh.model.geo.addLine(30, 31, 29)
    gmsh.model.geo.addLine(31, 32, 30)
    gmsh.model.geo.addLine(32, 29, 31)
    gmsh.model.geo.addLine(29, 30, 32)
    gmsh.model.geo.addCurveLoop([-29, -30, -31, -32], 108)

    gmsh.model.geo.addLine(34, 35, 33)
    gmsh.model.geo.addLine(35, 36, 34)
    gmsh.model.geo.addLine(36, 33, 35)
    gmsh.model.geo.addLine(33, 34, 36)
    gmsh.model.geo.addCurveLoop([-33, -34, -35, -36], 109)

    gmsh.model.geo.addLine(38, 39, 37)
    gmsh.model.geo.addLine(39, 40, 38)
    gmsh.model.geo.addLine(40, 37, 39)
    gmsh.model.geo.addLine(37, 38, 40)
    gmsh.model.geo.addCurveLoop([-37, -38, -39, -40], 110)

    gmsh.model.geo.addPlaneSurface([106, 107, 108, 109, 110], 2)

    # sides
    gmsh.model.geo.addLine(2, 22, 41)
    gmsh.model.geo.addLine(3, 23, 42)
    gmsh.model.geo.addLine(4, 24, 43)
    gmsh.model.geo.addLine(1, 21, 44)
    gmsh.model.geo.addLine(6, 26, 45)
    gmsh.model.geo.addLine(7, 27, 46)
    gmsh.model.geo.addLine(8, 28, 47)
    gmsh.model.geo.addLine(5, 25, 48)
    gmsh.model.geo.addLine(10, 30, 49)
    gmsh.model.geo.addLine(11, 31, 50)
    gmsh.model.geo.addLine(12, 32, 51)
    gmsh.model.geo.addLine(9, 29, 52)
    gmsh.model.geo.addLine(14, 34, 53)
    gmsh.model.geo.addLine(15, 35, 54)
    gmsh.model.geo.addLine(16, 36, 55)
    gmsh.model.geo.addLine(13, 33, 56)
    gmsh.model.geo.addLine(18, 38, 57)
    gmsh.model.geo.addLine(19, 39, 58)
    gmsh.model.geo.addLine(20, 40, 59)
    gmsh.model.geo.addLine(17, 37, 60)

    gmsh.model.geo.addCurveLoop([1, 42, -21, -41], 111)
    gmsh.model.geo.addPlaneSurface([-111], 3)

    gmsh.model.geo.addCurveLoop([2, 43, -22, -42], 112)
    gmsh.model.geo.addPlaneSurface([-112], 4)

    gmsh.model.geo.addCurveLoop([3, 44, -23, -43], 113)
    gmsh.model.geo.addPlaneSurface([-113], 5)

    gmsh.model.geo.addCurveLoop([4, 41, -24, -44], 114)
    gmsh.model.geo.addPlaneSurface([-114], 6)


    gmsh.model.geo.addCurveLoop([5, 46, -25, -45], 115)
    gmsh.model.geo.addPlaneSurface([-115], 7)

    gmsh.model.geo.addCurveLoop([6, 47, -26, -46], 116)
    gmsh.model.geo.addPlaneSurface([-116], 8)

    gmsh.model.geo.addCurveLoop([7, 48, -27, -47], 117)
    gmsh.model.geo.addPlaneSurface([-117], 9)

    gmsh.model.geo.addCurveLoop([8, 45, -28, -48], 118)
    gmsh.model.geo.addPlaneSurface([-118], 10)
    

    gmsh.model.geo.addCurveLoop([9, 50, -29, -49], 119)
    gmsh.model.geo.addPlaneSurface([-119], 11)

    gmsh.model.geo.addCurveLoop([10, 51, -30, -50], 120)
    gmsh.model.geo.addPlaneSurface([-120], 12)

    gmsh.model.geo.addCurveLoop([11, 52, -31, -51], 121)
    gmsh.model.geo.addPlaneSurface([-121], 13)

    gmsh.model.geo.addCurveLoop([12, 49, -32, -52], 122)
    gmsh.model.geo.addPlaneSurface([-122], 14)

    
    gmsh.model.geo.addCurveLoop([13, 54, -33, -53], 123)
    gmsh.model.geo.addPlaneSurface([-123], 15)

    gmsh.model.geo.addCurveLoop([14, 55, -34, -54], 124)
    gmsh.model.geo.addPlaneSurface([-124], 16)

    gmsh.model.geo.addCurveLoop([15, 56, -35, -55], 125)
    gmsh.model.geo.addPlaneSurface([-125], 17)

    gmsh.model.geo.addCurveLoop([16, 53, -36, -56], 126)
    gmsh.model.geo.addPlaneSurface([-126], 18)

    
    gmsh.model.geo.addCurveLoop([17, 58, -37, -57], 127)
    gmsh.model.geo.addPlaneSurface([127], 19)

    gmsh.model.geo.addCurveLoop([18, 59, -38, -58], 128)
    gmsh.model.geo.addPlaneSurface([128], 20)

    gmsh.model.geo.addCurveLoop([19, 60, -39, -59], 129)
    gmsh.model.geo.addPlaneSurface([129], 21)

    gmsh.model.geo.addCurveLoop([20, 57, -40, -60], 130)
    gmsh.model.geo.addPlaneSurface([130], 22)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end



"""
    meshcone(radius, height, h)

Create a mesh of an ice-cream cone of `radius` and `height`. The axis of symmetry is ̂z.

The target edge size is `h`.
"""
function meshcone(radius, height, h)
    fno = tempname() * ".msh"

    gmsh.initialize()
    gmsh.model.add("cone")
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 1)
    gmsh.model.geo.addPoint(0.0, 0.0, radius, h, 2)
    gmsh.model.geo.addPoint(0.0, 0.0, -height, h, 3)
    gmsh.model.geo.addPoint(radius, 0.0, 0.0, h, 4)
    gmsh.model.geo.addPoint(0.0, radius, 0.0, h, 5)
    gmsh.model.geo.addPoint(-radius, 0.0, 0.0, h, 6)
    gmsh.model.geo.addPoint(0.0, -radius, 0.0, h, 7)

    gmsh.model.geo.addCircleArc(4, 1, 5, 1)
    gmsh.model.geo.addCircleArc(5, 1, 6, 2)
    gmsh.model.geo.addCircleArc(6, 1, 7, 3)
    gmsh.model.geo.addCircleArc(7, 1, 4, 4)

    gmsh.model.geo.addCircleArc(2, 1, 5, 5)
    gmsh.model.geo.addCircleArc(2, 1, 6, 6)
    gmsh.model.geo.addCircleArc(2, 1, 7, 7)
    gmsh.model.geo.addCircleArc(2, 1, 4, 8)

    gmsh.model.geo.addLine(3, 4, 9)
    gmsh.model.geo.addLine(3, 5, 10)
    gmsh.model.geo.addLine(3, 6, 11)
    gmsh.model.geo.addLine(3, 7, 12)

    gmsh.model.geo.addCurveLoop([1, -5, 8], 1)
    gmsh.model.geo.addSurfaceFilling([1], 1)
    gmsh.model.geo.addCurveLoop([2, -6, 5], 2)
    gmsh.model.geo.addSurfaceFilling([2], 2)
    gmsh.model.geo.addCurveLoop([3, -7, 6], 3)
    gmsh.model.geo.addSurfaceFilling([3], 3)
    gmsh.model.geo.addCurveLoop([4, -8, 7], 4)
    gmsh.model.geo.addSurfaceFilling([4], 4)

    gmsh.model.geo.addCurveLoop([-1, 10, -9], 5)
    gmsh.model.geo.addSurfaceFilling([5], 5)
    gmsh.model.geo.addCurveLoop([-2, 11, -10], 6)
    gmsh.model.geo.addSurfaceFilling([6], 6)
    gmsh.model.geo.addCurveLoop([-3, 12, -11], 7)
    gmsh.model.geo.addSurfaceFilling([7], 7)
    gmsh.model.geo.addCurveLoop([-4, 9, -12], 8)
    gmsh.model.geo.addSurfaceFilling([8], 8)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end



"""
    meshsqpyramid(width, height, h)

Create a mesh of a square-based pyramid of `width` and `height`. The axis of symmetry is ̂z.

The target edge size is `h`.
"""
function meshsqpyramid(width, height, h)
    fno = tempname() * ".msh"

    gmsh.initialize()
    gmsh.model.add("sqpyramid")
    gmsh.model.geo.addPoint(0.0, 0.0, height, h, 1)
    gmsh.model.geo.addPoint(width/2, -width/2, 0.0, h, 2)
    gmsh.model.geo.addPoint(width/2, width/2, 0.0, h, 3)
    gmsh.model.geo.addPoint(-width/2, width/2, 0.0, h, 4)
    gmsh.model.geo.addPoint(-width/2, -width/2, 0.0, h, 5)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(1, 3, 2)
    gmsh.model.geo.addLine(1, 4, 3)
    gmsh.model.geo.addLine(1, 5, 4)

    gmsh.model.geo.addLine(2, 3, 5)
    gmsh.model.geo.addLine(3, 4, 6)
    gmsh.model.geo.addLine(4, 5, 7)
    gmsh.model.geo.addLine(5, 2, 8)

    gmsh.model.geo.addCurveLoop([-5, -6, -7, -8], 1)    
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop([1, -2, 5], 2)    
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addCurveLoop([2, -3, 6], 3)    
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addCurveLoop([3, -4, 7], 4)    
    gmsh.model.geo.addPlaneSurface([4], 4)
    gmsh.model.geo.addCurveLoop([4, -1, 8], 5)    
    gmsh.model.geo.addPlaneSurface([5], 5)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end



"""
    meshstarpyramid(majorradius, minorradius, nbofpoints, height, h)

Create a mesh of a star-based pyramid of star consisting of `nbofpoints` lying on `majorradius` and `nbofpoints` lying on `minorradius`, and `height`. The axis of symmetry is ̂z.

The target edge size is `h`.
"""
function meshstarpyramid(majorradius, minorradius, nbofpoints, height, h)
    fno = tempname() * ".msh"

    gmsh.initialize()
    gmsh.model.add("starpyramid")

    gmsh.model.geo.addPoint(0.0, 0.0, height, h, 0)

    angle = 2π/nbofpoints
    for i in 1:nbofpoints
        gmsh.model.geo.addPoint(majorradius*cos((i-1)*angle), majorradius*sin((i-1)*angle), 0.0, h, i)
        gmsh.model.geo.addPoint(minorradius*cos((i-0.5)*angle), minorradius*sin((i-0.5)*angle), 0.0, h, nbofpoints + i)
    end

    for i in 1:nbofpoints-1
        gmsh.model.geo.addLine(i, nbofpoints + i, i)
        gmsh.model.geo.addLine(nbofpoints + i, i + 1, nbofpoints + i)
    end

    # for case i = nbofpoints
    gmsh.model.geo.addLine(nbofpoints, 2*nbofpoints, nbofpoints)
    gmsh.model.geo.addLine(2*nbofpoints, 1, 2*nbofpoints)

    for i in 1:2*nbofpoints
        gmsh.model.geo.addLine(0, i, 2*nbofpoints + i)
    end

    base_curve_loop = Vector{Int}()
    for i in 1:nbofpoints
        append!(base_curve_loop, [-i, -(nbofpoints + i)])
    end

    @show base_curve_loop
    gmsh.model.geo.addCurveLoop(base_curve_loop, 3*nbofpoints)
    gmsh.model.geo.addPlaneSurface([3*nbofpoints], 3*nbofpoints)

    for i in 1:nbofpoints-1
        gmsh.model.geo.addCurveLoop([i, -(3*nbofpoints + i), 2*nbofpoints + i], i)
        gmsh.model.geo.addCurveLoop([nbofpoints + i, -(2*nbofpoints + i + 1), 3*nbofpoints + i], nbofpoints + i)
    end

    # for case i = nbofpoints
    gmsh.model.geo.addCurveLoop([nbofpoints, -4*nbofpoints, 3*nbofpoints], nbofpoints)
    gmsh.model.geo.addCurveLoop([2*nbofpoints, -(2*nbofpoints + 1), 4*nbofpoints], 2*nbofpoints)

    for i in 1:2*nbofpoints
        gmsh.model.geo.addPlaneSurface([i], i)
    end

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end