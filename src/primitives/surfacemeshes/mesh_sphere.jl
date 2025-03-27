using StaticArrays
# using Delaunay
using GmshTools

"""
    meshsphere(radius::F, h::F) where F

returns Mesh(vertices, faces)

Function gives a simplicial areal mesh of a sphere of radius and edge length.
Additionally, the function takes kwargs - generator, delaunay
-generator
:compsciencemeshes - is default
    -delaunay --- only a kwarg argument for compsciencemeshes
    :(2D) - is default, the triangulation is a 2D delaunay and is significantly faster
    than its 3D counterpart for medium to large meshes. The triangulation is done for
    x- and y- vertices, and the z- coordinate of the vertices are calculated from x, y 
    using Stereographic projection.
    :(3D) - calls the 3D Delaunay for triangulation, is significantly slower for 
    medium to large meshes.
:gmsh - returns a mesh using gmsh

The delaunay kwarg is only an argument for compsciencemeshes. The generator kwarg 
has precedence over delaunay.

Also see function - gmshsphere.
"""
function meshsphere(radius::F, h::F; 
    delaunay =:(2D), generator=:gmsh) where F
    if generator == :compsciencemeshes
        error("delaunay generator disabled waiting fix on Linux...")
        msh = mesh_sphere(radius, h; delaunay)
    elseif generator == :gmsh
        msh = gmshsphere(radius, h)
    else
        @error "generators are gmsh and compsciencemeshes only"
    end
    return msh
end

function meshsphere(;radius, h, delaunay=:(2D), generator=:gmsh)
    meshsphere(radius, h; delaunay=delaunay, generator=generator)
end

# Function to generate a uniform spherical mesh
function mesh_sphere(radius::F, len::F; delaunay =:(2D)) where F
    if delaunay ==:(3D)
        verts, faces = unit_centered_sphere(len/radius)
    elseif delaunay ==:(2D)
        verts, faces = unit_centered_sphere2(len/radius)
    end

    verts .= radius.*verts
    faces[(size(faces, 1) ÷ 2 + 1):end][2:3] .= 
        faces[(size(faces, 1) ÷ 2 + 1):end][3:-1:2]

    return Mesh(verts, faces)
end

# Function to generate a unit sphere mesh using delaunay 3D function
function unit_centered_sphere(dist::F) where F
    #steps down from the top to list vertices of the hemisphere
    thmax = F(π/2)
    points = []
    #number of triangles the azimuthal/ϕ can accomodate
    nth = Int(floor(thmax / F(dist * sqrt(3) / 2))) + 1
    th_values = collect(LinRange(F(0), thmax, nth))
    upperpts = 0
    for th in th_values
        if th == F(π/2)
            upperpts = size(points, 1)
        end
        r = sin(th)
        #number of triangles the circle at ϕ can accomodate
        nph = Int(floor(2 * π * r / dist)) + 1
        #eps() to avoid error in alignment
        ph = collect(LinRange(2*π/nph, 2*π, nph)) 
            .+ (π / (th + eps(F)))^2 
        for (i, θ) in enumerate(ph)
            x = r * cos(θ)
            y = r * sin(θ)
            #vertices (x, y, z)
            append!(points, [
                [x, 
                y, 
                sign(pi/2-th)*F(sqrt(1-min(x.^2+y.^2,1)))]
                ])
        end       
    end

    # introducing a dummy point
    append!(points, [[F(0), F(0), F(0)]])
    #total number of vertices (x, y, z) in the sphere
    V = length(points) + upperpts 
    sp = reshape(
        [point[i] for i in 1:3 for point in points], 
        (length(points), 3)
        )
    verts = zeros(SVector{3, F}, V - 1)
    #delaunay 3D returns tetrahedrons with the fourth point as the dummy point
    tris = delaunay_triangulation(sp[:, 1:3]).simplices
    #removing the dummy point
    pop!(points)

    #assigning vertices
    for i in eachindex(points)
        verts[i] = [points[i][1], points[i][2], points[i][3]]
        if i <= upperpts
            verts[length(points) + i] = 
            [verts[i][1], verts[i][2], -verts[i][3]]
        end
    end

    #for the bottom hemisphere
    vertreindex = collect(1:V - 1) 
    vertreindex[1:upperpts] = V .- reverse(vertreindex[1:upperpts])
    temp = []
    #to remove the fourth vertex/dummy point from all faces
    for i in 1:size(tris, 1)
        id = findall(x -> x == V - upperpts, tris[i, :])
        if id != []
            append!(temp, [deleteat!(tris[i, :], id[1])])
        end
    end
    for i in 1:length(temp)
        append!(temp, [vertreindex[temp[i]]])
    end
    faces = zeros(SVector{3, Int}, length(temp))
    #faces 
    for i in 1:length(temp)
       faces[i] = SVector{3, Int}(temp[i])
    end
    return verts, faces
end


#Function to generate a unit sphere mesh using delaunay 2D function
function unit_centered_sphere2(dist::F) where F
    #steps down from the top to list vertices of the hemisphere
    thmax = F(π / 2)
    points = []
    nth = Int(floor(thmax / F(dist * sqrt(3) / 2))) + 1
    th_values = collect(LinRange(F(0), thmax, nth))
    upperpts = 0
    for th in th_values
        if th == F(π / 2)
            upperpts = size(points, 1)
        end
        r = sin(th)
        nph = Int(floor(2 * π * r / dist)) + 1
        ph = collect(LinRange(2*π/nph, 2*π, nph)) 
            .+ (π / (th + eps(F)))^2
        for (i, θ) in enumerate(ph)
            x = th * cos(θ)
            y = th * sin(θ)
            #(x, y) coordinates of the vertices
            append!(points, [[x, y]])
        end       
    end
    #total number of vertices (x, y, z)
    V = length(points) + upperpts 
    sp = reshape(
        [point[i] for i in 1:2 for point in points], 
        (length(points), 2)
        )
    verts = zeros(SVector{3, F}, V)
    #finds the norm of the coordinates
    t = leonorm(points)
    mfact = sinc.(t)
    #delaunay 2D returns triangles
    tris = delaunay_triangulation(sp[:, 1:2]).simplices
    #vertices (x, y, z)
    for i in eachindex(points)
        verts[i] = [
            points[i][1]*mfact[i], 
            points[i][2]*mfact[i], 
            cos(t[i])
            ]
        if i <= upperpts
            verts[length(points) + i] = [
                verts[i][1], verts[i][2], -verts[i][3]
                ]
        end
    end
    #to find the vertices on the bottom hemisphere
    vertreindex = collect(1:V)
    vertreindex[1:upperpts] = V + 1 .- reverse(
        vertreindex[1:upperpts]
        )
    tris = vcat(tris, vertreindex[tris])
    #faces
    faces = zeros(SVector{3, Int}, size(tris, 1))
    for i in 1:size(tris, 1)
        faces[i] = SVector{3, Int}(tris[i, 1:3])
    end
    return verts, faces
end

#Function to compute norm - used by unitCenteredSphere2
function leonorm(inmat)
    outmat = zeros(length(inmat))
    #a column vector of summed squares of (x,y)
    for i in 1:length(inmat)
        #x^2 + y^2 along each row
        outmat[i] = inmat[i][1]*conj(inmat[i][1]) + inmat[i][2]*conj(inmat[i][2])
    end
    return sqrt.(outmat)
end

# Function to compute sin(x)/x - used by unitCenteredSphere2
function sinc(x::F) where F
    if x == F(0)
        result = F(1)
    else
        result = sin.(x) ./ x
    end
    return result
end


# Function to perform Delaunay triangulation using Delaunay package
function delaunay_triangulation(points)
    triangulation = delaunay(points)
    return triangulation  
end

"""
    gmshsphere(radius, delta)

Create a mesh of a sphere of radius `radius` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

The target edge size is `delta`.
"""
function gmshsphere(radius, delta; tempname=tempname())
    s = """
lc = $delta;

Point(1)={0,0,0,lc};
Point(2)={$radius,0,0,lc};
Point(3)={0,$radius,0,lc};
Point(4)={-$radius,0,0,lc};
Point(5)={0,-$radius,0,lc};
Point(6)={0,0,$radius,lc};
Point(7)={0,0,-$radius,lc};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};
Circle(5)={6,1,3};
Circle(6)={3,1,7};
Circle(7)={7,1,5};
Circle(8)={5,1,6};

Line Loop(1)={-1,-6,-7,-4};
Line Loop(2)={1,-5,-8,4};
Line Loop(3)={-2,-3,7,6};
Line Loop(4)={2,3,8,5};

Ruled Surface(1)={1} In Sphere{1};
Ruled Surface(2)={2} In Sphere{1};
Ruled Surface(3)={3} In Sphere{1};
Ruled Surface(4)={4} In Sphere{1};
"""

    fn = tempname
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end
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