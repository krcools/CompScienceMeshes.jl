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