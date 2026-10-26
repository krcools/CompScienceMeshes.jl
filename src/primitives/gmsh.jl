using GmshTools

"""
    meshgeo(geofile; physical=nothing, dim=2, tempname=tempname(), kwargs...)

"""
function meshgeo(geofile; physical=nothing, dim=2, order=1, tempname=tempname(), kwargs...)
    
    io = open(geofile)
    str = read(io, String)
    close(io)
    
    for (key,val) in kwargs
        key_str = String(key)
        pat = Regex("$key_str\\s*=\\s*[0-9]*.?[0-9]*;")
        sub = SubstitutionString("$key_str = $val;")
        str = replace(str, pat => sub)
    end

    # println(str)
    # return
    
    temp_geo = tempname
    open(temp_geo, "w") do io
        print(io, str)
    end

    temp_msh = tempname * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(temp_geo)
    gmsh.model.mesh.generate(dim)
    order > 1 && gmsh.model.mesh.setOrder(order)
    gmsh.write(temp_msh)
    gmsh.finalize()

    # @show temp_msh

    # run(`gmsh $temp_geo -2 -format msh2 -o $temp_msh`)
    if dim == 2
        # m = read_gmsh_mesh(temp_msh, physical=physical)
        m = load_gmsh_mesh(temp_msh; element=:triangle, physical, order)
    elseif dim == 3
        # m = read_gmsh3d_mesh(temp_msh, physical=physical)
        m = load_gmsh_mesh(temp_msh; element=:tetrahedron, physical, order)
    else
        error("gmsh files of dimension $(dim) cannot be read.")
    end
    
    rm(temp_msh)
    rm(temp_geo)
    
    return m
end