using GmshTools

"""
    meshdisk(radius::T, delta::T, tempname=tempname()) where T<:Real

"""
function meshdisk(radius::T, delta::T, tempname=tempname()) where T<:Real
    s = """
lc = $delta;

Point(1)={0,0,0,lc};
Point(2)={$radius,0,0,lc};
Point(3)={0,$radius,0,lc};
Point(4)={-$radius,0,0,lc};
Point(5)={0,-$radius,0,lc};


Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};


Line Loop(1)={2,3,4,1};

Plane Surface(1) = {1};
"""

    fn = tempname
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname * ".msh"

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(fn)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)

    rm(fno)
    rm(fn)

    return m

end
