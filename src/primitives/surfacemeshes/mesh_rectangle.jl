using StaticArrays
using GmshTools

"""
    meshrectangle(len::F, breadth::F, edge_length::F, udim = 3) where F

returns Mesh(vertices, faces)

Function calls arg - udim or universal dimension -
    2 - returns a rectangle in a 2D space
    3 - is default, returns a rectangle with the 3D coordinate system at z = 0

Function also calls a kwarg - generator - 
:compsciencemeshes - is default, is pre-allocated with its vectors 
    of nodes and faces. It is a structured mesh. 

    Nodes:
    (m+1) nodes along b
    (n+1) nodes along a
    = (m+1)*(n+1) total nodes

    Faces:
    m elements along b 
    n elements along a 
    = 2*m*n total faces -> 2 triangles in each rectangular face

    The faces along y-axis can be obtained by (2*element number) and 
    (2*element number - 1), this is repeated for each element in x-axis.

    The odd faces are the triangles right-angled at bottom,
      /|
     /_|      or laterally-inverted

     the even at top
      _      
     | /    
     |/        or laterally-inverted

:gmsh - returns an unstructured mesh using gmsh, which in turn using 
Delaunay triangulation. 

kwarg (not implemented yet): boundary_only - returns the mesh of the boundary 
of the rectangle, if true

Also see gmsh function - gmshrectangle.
"""
function meshrectangle(len::F, breadth::F, edge_length::F, udim = 3; 
    boundary_only = false, generator = :compsciencemeshes, element=:triangle) where F
    if generator == :gmsh && element == :triangle
        @assert udim==3 
        msh = gmshrectangle(len, breadth, edge_length)
    elseif generator == :compsciencemeshes && element == :quadrilateral
        @assert udim==3
        msh =  meshrectangle_quadrilateral(len, breadth, edge_length)
    elseif generator == :compsciencemeshes && element == :triangle
        @info "Generating a structured mesh: The dimensions of the rectangle are 
            approximated by multiples of edge length.
            For exact dimensions/ unstructured grid, use kwarg - generator = :gmsh"
        msh = mesh_rectangle(len, breadth, edge_length, udim)
    else 
        @error "generators are gmsh and compsciencemeshes only"
    end

    if boundary_only == true
        @error "not implemented yet"
    end
    return msh
end 

"""
    mesh_rectangle(a::Float64, b::Float64, h::Float64, udim)

    The mesh function is pre-allocated with its vectors of nodes and faces. It is a 
    structured mesh.
"""
function mesh_rectangle(a::F, b::F, h::F, udim) where F
    @assert udim == 2 || udim == 3 "Universal dimension can only be 2 or 3"
    #structured mesh:  isapprox(a%h, F(0)) && isapprox(b%h, F(0))
    n = Int(round(a/h))  # number of elements along a
    m = Int(round(b/h))  # number of elements along b
    
    nodes = zeros(SVector{udim, F}, (m + 1)*(n + 1))
    faces = Vector{SVector{3, Int64}}(undef, 2*m*n)
    
    for ix in range(0, n - 1)
        for iy in range(1, m)
            if udim == 3
                node = SVector((ix)*h, (iy - 1)*h, F(0))
            else
                node = SVector((ix)*h, (iy - 1)*h)
            end
            
            nodes[(ix)*(m + 1) + iy] = node
            face = SVector(
                (ix)*(m + 1) + (iy),
                (ix)*(m + 1) + (iy + 1),
                (ix + 1)*(m + 1) + (iy)
            )
            faces[(ix)*2*m + (2*iy - 1)] = face
            
            face = SVector(
                (ix)*(m + 1) + (iy + 1), 
                (ix + 1)*(m + 1) + (iy + 1), 
                (ix + 1)*(m + 1) + (iy)
                )
            faces[(ix)*2*m + (2*iy)] = face
        end
        # for the mth element in y-direction
        if udim == 3
            nodes[(ix + 1)*(m + 1)] = SVector((ix)*h, m*h, F(0))
        else
            nodes[(ix + 1)*(m + 1)] = SVector((ix)*h, m*h)
        end
        
    end
    # for ix = n
    for iy in range(0, m)
        if udim == 3
            nodes[n*(m + 1) + iy + 1] = SVector(n*h, (iy*h), F(0))
        else 
            nodes[n*(m + 1) + iy + 1] = SVector(n*h, (iy*h))
        end
    end
        
    return Mesh(nodes, faces)
end

"""
    gmshrectangle(width, height, delta)
	Meshes unstructured rectangle (Delaunay Triangulation).
    Takes kwarg - tempname for naming the .geo file
"""
function gmshrectangle(width, height, delta; tempname=tempname())
    s =
		"""
		lc = $delta;

		Point(1)={0,0,0,lc};
		Point(4)={$width,0,0,lc};
		Point(5)={0,$height,0,lc};
		Point(8)={$width,$height,0,lc};

		Line(4)={4,1};
		Line(8)={8,5};
		Line(9)={1,5};
		Line(12)={4,8};

		Line Loop(5)={4,-12,-8,9};

		Plane Surface(5)={5};
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
    m = meshrectangle(2.0, 2.0, 1.0; element=:quadrilateral)
    @test length(m) == 4
    ch = chart(m, first(m))
    p = neighborhood(ch, point(0.5, 0.5))
    @test cartesian(p) ≈ point(0.5, 0.5, 0.0)
end
