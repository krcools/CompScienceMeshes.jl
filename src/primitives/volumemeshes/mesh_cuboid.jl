using StaticArrays

"""
    tetmeshcuboid(length::F, breadth::F, width::F, edge length::F) where F

returns Mesh(vertices, faces)

Function returns a simplicial volumetric mesh of a cuboid. It takes kwarg - 
generator 
:compsciencemeshes - is default, and gives a structured tetrahedral mesh of 
    a cuboid, the dimensions of the cuboid are approximated by multiples 
    of edge length -

    there are 24 tetrahedrons in each unit cube - 4 on each facet - 
    and there are n*m*p unit cubes in each cuboid, where n, m, p are number of
    segments along length, breadth, width respectively.
    
:gmsh - returns a tetrahedral mesh of a cuboid using gmsh. 

Also see the gmsh function - tetgmshcuboid.
"""
function tetmeshcuboid(len::F, breadth::F, width::F, edge_len::F; 
    generator = :compsciencemeshes) where F
    if generator == :compsciencemeshes
        # @info "Generating a structured mesh: The dimensions of the cuboid are 
        #     approximated by multiples of edge length.
        #     For exact dimensions/ unstructured grid, use kwarg - generator = :gmsh"
        msh = tetmesh_cuboid(len, breadth, width, edge_len)
    elseif generator == :gmsh
        msh = tetgmshcuboid(len, breadth, width, edge_len)
    else
        @error "generators are gmsh and compsciencemeshes only"
    end
    return msh
end

@generated function tetmesh_cuboid(a, b, c, h)
    Core.println("Generating a structured mesh: The dimensions of the cuboid are 
            approximated by multiples of edge length. For exact dimensions and
            unstructured grids, use kwarg - generator = :gmsh")
    return :(tetmesh_cuboid_impl(a, b, c, h))
end

function tetmesh_cuboid_impl(a::F, b::F, c::F, h::F) where F
    n = Int(round(a/h)) #number of elements along x
    m = Int(round(b/h)) #number of elements along y
    p = Int(round(c/h)) #number of elements along z

    #total number of vertices becomes (m + 1) along each edge in y direction 
        #and along each plane in z direction
    vertices = Vector{SVector{3, F}}(
        undef, 
        (m + 1)*(n + 1)*(p + 1) + 4*m*n*p + (m*n + n*p + m*p)
        )
    #there are 24 tetrahedrons in each unit cube - 4 on each facet - 
        #and there are m*n*p unit cubes in each cuboid
    faces = Vector{SVector{4, Int}}(undef, 24*m*n*p)

    verts = (m + 1)*(n + 1)*(p + 1)
    cents = m*n*p
    sides = verts + cents

    for iz in 1 : p+1
        for iy in 1 : m+1
            for ix in 1 : n+1
                v = Vector{SVector{4, Int}}(undef, 24)
                #vertices of unit cubes
                ver_mem1 = ((ix - 1)*(m + 1) 
                         + (iy) 
                         + (iz - 1)*(m + 1)*(n + 1))
                vertices[ver_mem1] = SVector(
                    (ix - 1)*h, 
                    (iy - 1)*h, 
                    (iz - 1)*h
                    ) 
                #centers on the surfaces of unit cubes
                #they are ordered from front to back, bottom to top, left to right 
                if (ix != n + 1)&&(iy != m + 1)
                    ver_mem21 = (sides 
                              + iy 
                              + (ix - 1)*m 
                              + (iz - 1)*m*n)
                    vertices[ver_mem21] = SVector(
                        (ix - 0.5)*h, 
                        (iy - 0.5)*h, 
                        (iz - 1)*h
                        )
                end
                if (ix != (n + 1))&&(iz != (p + 1))
                    ver_mem22 = (sides + m*n*(p + 1) 
                              + iz 
                              + (ix - 1)*p 
                              + (iy - 1)*n*p)
                    vertices[ver_mem22] = SVector(
                        (ix - 0.5)*h, 
                        (iy - 1)*h, 
                        (iz - 0.5)*h
                        )
                end
                if (iy != (m + 1))&&(iz != (p + 1))
                    ver_mem23 = (sides + n*p*(m + 1) + m*n*(p + 1) 
                              + iy 
                              + (iz - 1)*m 
                              + (ix - 1)*m*p)
                    vertices[ver_mem23] = SVector(
                        (ix - 1)*h, 
                        (iy - 0.5)*h, 
                        (iz - 0.5)*h
                        )
                end

                if (ix != (n + 1))&&(iy != (m + 1))&&(iz != (p + 1))
                    #centers of unit cubes
                    ver_mem3 = (verts 
                             + (ix - 1)*m 
                             + iy 
                             + (iz - 1)*m*n)
                    vertices[ver_mem3] = SVector(
                        (ix - 0.5)*h, 
                        (iy - 0.5)*h, 
                        (iz - 0.5)*h
                        )

                    #faces on each unit cube
                    face_mem = ((ix - 1)*m + iy + (iz - 1)*m*n) - 1
                    #faces
                    for l in 0:2
                        for k in 0:1
                            for i in 0:1
                                for j in 0:1
                                    v1 = ((ix - 1 + i)*(m + 1) 
                                       + (iy + j) 
                                       + (iz - 1 + k)*(m + 1)*(n + 1))
                                    if l == 0 #front-back
                                        v2 = ((ix - 1 + j)*(m + 1) 
                                           + (iy + (1 - i)) 
                                           + (iz - 1 + k)*(m + 1)*(n + 1))
                                        v3 = (sides 
                                           + iy 
                                           + (ix - 1)*m 
                                           + (iz - 1 + k)*m*n)
                                    elseif l == 1 #bottom-top
                                        v2 = ((ix - 1 + k)*(m + 1) 
                                           + (iy + j) 
                                           + (iz - 1 + (1 - i))*(m + 1)*(n + 1))
                                        v3 = (sides + m*n*(p + 1) 
                                           + iz 
                                           + (ix - 1)*p 
                                           + (iy - 1 + j)*n*p)
                                    else #left-right
                                        v2 = ((ix - 1 + i)*(m + 1) 
                                           + (iy + (1 - k)) 
                                           + (iz - 1 + j)*(m + 1)*(n + 1))
                                        v3 = (sides + n*p*(m + 1) + m*n*(p + 1) 
                                           + iy + (iz - 1)*m 
                                           + (ix - 1 + i)*m*p)
                                    end
                                    v4 = ver_mem3 #always the center of the unit cube
                                    v[8*l + 4*k + 2*i + j + 1] = SVector(v1, v2, v3, v4)
                                end
                            end
                        end
                    end
                    for i in 1:24
                        faces[face_mem*24 + i] = v[i]
                    end

                end
            end
        end
    end

    return Mesh(vertices, faces)
end

"""
    tetgmshcuboid(width, height, length, delta)

returns a mesh of tetrahedral mesh of a cuboid using gmsh. The name of the 
.geo file can be specified with the kwarg - tempname.
"""
function tetgmshcuboid(width, height, length, delta; tempname = tempname())
    s =
"""
lc = $delta;

Point(1)={0,0,0,lc};
Point(2)={0,0,$length,lc};
Point(3)={$width,0,$length,lc};
Point(4)={$width,0,0,lc};
Point(5)={0,$height,0,lc};
Point(6)={0,$height,$length,lc};
Point(7)={$width,$height,$length,lc};
Point(8)={$width,$height,0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,5};
Line(9)={1,5};
Line(10)={2,6};
Line(11)={3,7};
Line(12)={4,8};

Line Loop(1)={-1,-2,-3,-4};
Line Loop(2)={1,-9,-5,10};
Line Loop(3)={2,-10,-6,11};
Line Loop(4)={3,-11,-7,12};
Line Loop(5)={4,-12,-8,9};
Line Loop(6)={5,6,7,8};

Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};
Plane Surface(4)={4};
Plane Surface(5)={5};
Plane Surface(6)={6};

Surface Loop(1)={1,2,3,4,5,6};

Volume(1)={1};
Physical Volume(1) = {1};
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
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(fn)
    gmsh.model.mesh.generate(3)
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh3d_mesh(fno)

    rm(fno)
    rm(fn)
    return m

end