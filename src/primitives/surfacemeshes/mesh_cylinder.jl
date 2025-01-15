using GmshTools

"""
    meshcylinder(;radius, height, h, tempname=tempname())

"""
function meshcylinder(;radius, height, h, tempname=tempname())

    fn = joinpath(@__DIR__,"geos/cylinder3.geo")
    io = open(fn)
    str = read(io, String)
    close(io)

    str = replace(str, "r = 1.0;" => "r = $radius;")
    str = replace(str, "z = 1.0;" => "z = $height;")
    str = replace(str, "h = 1.0;" => "h = $h;")

    temp_geo = tempname
    open(temp_geo, "w") do io
        print(io, str)
    end

    temp_msh = tempname * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(temp_geo)
    gmsh.model.mesh.generate(3)
    gmsh.write(temp_msh)
    gmsh.finalize()
    m = read_gmsh3d_mesh(temp_msh)

    rm(temp_msh)
    rm(temp_geo)

    return m
end

