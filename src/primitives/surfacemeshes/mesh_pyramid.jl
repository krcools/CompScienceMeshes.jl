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