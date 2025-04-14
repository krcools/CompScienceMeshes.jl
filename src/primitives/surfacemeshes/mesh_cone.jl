"""
    meshcone(radius, height, h)

Create a mesh of an ice-cream cone of `radius` and `height`. The axis of symmetry is ̂z.

The target edge size is `h`.
"""
function meshcone(radius, height, h)
    fno = tempname() * ".msh"

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
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