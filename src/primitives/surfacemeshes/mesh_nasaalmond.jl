

"""
    meshnasaalmond(h)

Create a mesh of the NASA almond - the RCS benchmark. 
The NASA almond geometry is described in: Woo et. al., IEEE Antennas and Propagation Magazine, Vol. 35, No. 1, p. 84-89, 1993. 

for t1 < t < 0, -pi < psi < pi
    x = d*t
    y = y1*d*sqrt(1 - (t/t1)^2)*cos(psi)
    z = z1*d*sqrt(1 - (t/t1)^2)*sin(psi)

for 0 < t < t2, -pi < psi < pi
    x = d*t
    y = y2*d*(sqrt(1 - (t/t2)^2) - 0.96)*cos(psi)
    z = z2*d*(sqrt(1 - (t/t2)^2) - 0.96)*sin(psi)

The target edge size is `h`.
"""

function meshnasaalmond(h)

    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("nasaalmond")

    # parameters
    d = 9.936
    t1= -0.416667
    t2 = 1 + t1
    x1, x2 = -t1, t2/0.28
    y1, z1 = 0.193333, 0.064444
    y2, z2 = 25*y1, 25*z1

    fxy1 = -sqrt(x1^2 - y1^2)
    fxz1 = -sqrt(x1^2 - z1^2)
    fyz1 = sqrt(y1^2 - z1^2)

    # First part of the surface
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 0)
    
    gmsh.model.geo.addPoint(t1*d, 0.0, 0.0, h, 1)
    
    gmsh.model.geo.addPoint(0.0, y1*d, 0.0, h, 2)
    gmsh.model.geo.addPoint(0.0, -y1*d, 0.0, h, 3)
    gmsh.model.geo.addPoint(fxy1, 0.0, 0.0, h, 4)

    gmsh.model.geo.addEllipseArc(3, 0, 4, 1, 1)
    gmsh.model.geo.addEllipseArc(1, 0, 4, 2, 2)
    

    gmsh.model.geo.addPoint(0.0, 0.0, z1*d, h, 5)
    gmsh.model.geo.addPoint(0.0, 0.0, -z1*d, h, 6)
    gmsh.model.geo.addPoint(fxz1, 0.0, 0.0, h, 7)

    gmsh.model.geo.addEllipseArc(6, 0, 7, 1, 3)
    gmsh.model.geo.addEllipseArc(1, 0, 7, 5, 4)

    gmsh.model.geo.addPoint(0.0, fyz1, 0.0, h, 8)
    gmsh.model.geo.addPoint(0.0, -fyz1, 0.0, h, 9)

    gmsh.model.geo.addEllipseArc(5, 0, 8, 2, 5)
    gmsh.model.geo.addEllipseArc(2, 0, 8, 6, 6)
    gmsh.model.geo.addEllipseArc(6, 0, 9, 3, 7)
    gmsh.model.geo.addEllipseArc(3, 0, 9, 5, 8)

    gmsh.model.geo.addCurveLoop([-1, -4, 8], 1)
    gmsh.model.geo.addSurfaceFilling([1], 1)
    gmsh.model.geo.addCurveLoop([4, -2, 5], 2)
    gmsh.model.geo.addSurfaceFilling([2], 2)
    gmsh.model.geo.addCurveLoop([2, 3, 6], 3)
    gmsh.model.geo.addSurfaceFilling([3], 3)
    gmsh.model.geo.addCurveLoop([-3, 1, 7], 4)
    gmsh.model.geo.addSurfaceFilling([4], 4)

    # Second part of the surface
    gmsh.model.geo.addPoint(t2*d, 0.0, 0.0, h, 10)

    Np = 10                                                 # number of sampling points on the curves
    yc1, yc2, zc1, zc2 = [10], [10], [10], [10]             # points for the curves

    for i = 1:Np-1
        ti = t2/Np*(Np-i)
        sc = sqrt(1 - (ti/x2)^2) - 0.96

        gmsh.model.geo.addPoint(ti*d, y2*d*sc, 0.0, h, 10*i+1)
        gmsh.model.geo.addPoint(ti*d, -y2*d*sc, 0.0, h, 10*i+2)
        gmsh.model.geo.addPoint(ti*d, 0.0, z2*d*sc, h, 10*i+3)
        gmsh.model.geo.addPoint(ti*d, 0.0, -z2*d*sc, h, 10*i+4)

        append!(yc1, [10*i+1])
        append!(yc2, [10*i+2])
        append!(zc1, [10*i+3])
        append!(zc2, [10*i+4])
    end

    append!(yc1, [2])
    append!(yc2, [3])
    append!(zc1, [5])
    append!(zc2, [6])

    gmsh.model.geo.addBSpline(yc1, 9)
    gmsh.model.geo.addBSpline(yc2, 10)
    gmsh.model.geo.addBSpline(zc1, 11)
    gmsh.model.geo.addBSpline(zc2, 12)

    gmsh.model.geo.addCurveLoop([-5, 9, -11], 5)
    gmsh.model.geo.addSurfaceFilling([5], 5)
    gmsh.model.geo.addCurveLoop([-6, -9, 12], 6)
    gmsh.model.geo.addSurfaceFilling([6], 6)
    gmsh.model.geo.addCurveLoop([-7, -12, 10], 7)
    gmsh.model.geo.addSurfaceFilling([7], 7)
    gmsh.model.geo.addCurveLoop([-8, -10, 11], 8)
    gmsh.model.geo.addSurfaceFilling([8], 8)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    m = CompScienceMeshes.read_gmsh_mesh(fno)
    rm(fno)
    return m
end

