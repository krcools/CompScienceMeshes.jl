using StaticArrays
using GmshTools

"""
    meshcuboid(length::F, breadth::F, width::F, edge length::F) where F

returns Mesh(vertices, faces)

Function returns a simplicial areal mesh of a cuboid. It takes kwarg - generator 
:compsciencemeshes - is default, and gives a structured mesh of 
    a cuboid, the dimensions of the cuboid are approximated by multiples 
    of edge length -
    
    Number of faces = 2*(2*m*n) + 2*(2*n*p) + 2*(2*m*p) -> front, back; top, bottom;
    left, right
    where m is the number of elements along y, n along x, p along z

    The odd faces are the triangles right-angled at bottom,
      /|
     /_|      or laterally-inverted

     the even at top
      _      
     | /    
     |/        or laterally-inverted

    Number of vertices is 2*(n + 1)*(m + 1) + 2*(m + 1)*(p - 1) + 2*(n - 1)*(p - 1)
    to avoid repeating the number of nodes along the edges of the cuboid; 
    all nodes along front and back are stacked, 
    then the nodes along the left and right, skipping the nodes along the left-front, 
    left-back, right-front and right-back
    then the nodes along the top and bottom, skipping the nodes along the top-front, 
    -back, -left, -right, and bottom-front, -back, -left, -right.
    
    Since, it is a hollow cuboidal mesh,
    the nodes are numbered all along the front and back faces, first, and then, 
    along the x-y edges in z - direction, like so

                                 20----- 23----- 26
                                 |       |       |
                                 19      22      25
    back =>                      |       |       |
                                 18----- 21----- 24


                                 12----- 14----- 17
    middle                       |               |
    layers =>                    11      .       16
                                 |               |
                                 10----- 13----- 15


                                 3 ----- 6 ----- 9
                                 |       |       |
                                 2       5       8
    front =>                     |       |       |
                                 1 ----- 4 ----- 7

:gmsh - returns a areal mesh of a cuboid using gmsh and Closed Box only.

For other physical configurations of the surface of the cube, see gmshcuboid function.
"""
function meshcuboid(len::F, breadth::F, width::F, edge_len::F; 
    generator = :compsciencemeshes) where F
    if generator == :gmsh
        msh = gmshcuboid(len, breadth, width, edge_len)
    elseif  generator == :compsciencemeshes
        @info "Generating a structured mesh: The dimensions of the cuboid are 
            approximated by multiples of edge length.
            For exact dimensions/ unstructured grid, use kwarg - generator = :gmsh"
        msh = mesh_cuboid(len, width, width, edge_len)
    else
        @error "generators are gmsh and compsciencemeshes only"
    end
    return msh
end

#code for meshing a cuboid regularly
"""
    mesh_cuboid(a::F, b::F, c::F, h::F)

returns an areal structured mesh of a cuboid.

    """
function mesh_cuboid(a::F, b::F, c::F, h::F) where F
    # if  isapprox(a%h, F(0)) && isapprox(b%h, F(0) && isapprox(c%h, F(0)))
        n = Int(round(a/h))  # number of elements along a
        m = Int(round(b/h))  # number of elements along b
        p = Int(round(c/h))  #number of elements along c        

        nodes = Vector{SVector{3, F}}(undef, 2*(m*n + m*p + n*p + 1))
        faces = Vector{SVector{3, Int64}}(undef, 4*m*n + 4*n*p + 4*m*p)

        #along x-y 
        #the first face/node is the variable + 1
        back_node = (m + 1)*(n + 1) + 2*(p - 1)*(m + n)
        back_face = 2*m*p + 2*n*p + 2*m*n
        for ix in range(1, n)
            for iy in range(1, m)
                #nodes
                #front
                nodes[(ix - 1)*(m + 1) + iy] = SVector((ix - 1)*h, (iy - 1)*h, F(0))
                #back
                nodes[back_node + (ix - 1)*(m + 1) + iy] = SVector(
                    (ix - 1)*h, 
                    (iy - 1)*h, 
                    p*h
                    )
                #faces
                #front 
                faces[(ix - 1)*2*m + (2*iy - 1)] = SVector(
                    (ix - 1)*(m + 1) + (iy),
                    (ix)*(m + 1) + (iy),
                    (ix - 1)*(m + 1) + (iy + 1)
                )
                faces[(ix - 1)*2*m + (2*iy)] = SVector(
                    (ix - 1)*(m + 1) + (iy + 1), 
                    (ix)*(m + 1) + (iy), 
                    (ix)*(m + 1) + (iy + 1)
                )
                #back
                faces[back_face + (ix - 1)*2*m + (2*iy - 1)] = SVector(
                    back_node + (ix - 1)*(m + 1) + (iy), 
                    back_node + (ix - 1)*(m + 1) + (iy + 1), 
                    back_node + (ix)*(m + 1) + (iy)
                )
                faces[back_face + (ix - 1)*2*m + (2*iy)] = SVector(
                    back_node + (ix - 1)*(m + 1) + (iy + 1), 
                    back_node + (ix)*(m + 1) + (iy + 1), 
                    back_node + (ix)*(m + 1) + (iy)
                )                
            end
            # for the mth element in y-direction
            #front
            nodes[ix*(m + 1)] = SVector((ix - 1)*h, m*h,  F(0))
            #back
            nodes[back_node + ix*(m + 1)] = SVector(
                (ix - 1)*h, m*h, p*h)
        end
        # for ix = n
        for iy in range(0, m)
            #front
            nodes[n*(m + 1) + iy + 1] = SVector(n*h, (iy*h),  F(0))
            #back
            nodes[back_node + n*(m + 1) + iy + 1] = SVector(
                n*h, (iy*h), p*h) 
        end 

        # along y - z
        for iy in range(1, m + 1)
            for iz in range(2, p)
                #left node + m + 2*n - 1 gives the right nodes
                left_node_front = (n + 1)*(m + 1) + (iz - 2)*(2*(m + n))
                left_node_back = (n + 1)*(m + 1) + (iz - 3)*(2*(m + n))
                right_face = 2*m*n + 2*n*p + (iz - 2)*2*m
                left_face = 4*m*n + 4*n*p + 2*m*p + (iz - 2)*2*m
                #ix in 1 -> left nodes
                nodes[left_node_front + iy] = SVector(
                    F(0), (iy - 1)*h, (iz - 1)*h)
                if iz == 2 && iy != (m + 1)
                    #left faces
                    faces[left_face + (2*iy - 1)] = SVector(
                        iy, 
                        iy + 1, 
                        left_node_front + iy
                        )
                    faces[left_face + (2*iy)] = SVector(
                        iy + 1, 
                        left_node_front + iy + 1, 
                        left_node_front + iy
                        )
                    #right faces
                    faces[right_face + (2*iy - 1)] = SVector(
                        n*(m + 1) + iy, 
                        left_node_front + m + 2*n - 1 + iy, 
                        n*(m + 1) + iy + 1
                        )
                    faces[right_face + (2*iy)] = SVector(
                        n*(m + 1) + iy + 1, 
                        left_node_front + m + 2*n - 1 + iy, 
                        left_node_front + m + 2*n - 1 + iy + 1
                        )
                elseif iy != (m + 1)
                    #left faces
                    faces[left_face + (2*iy - 1)] = SVector(
                        left_node_back + iy, 
                        left_node_back + iy + 1, 
                        left_node_front + iy
                        )
                    faces[left_face + (2*iy)] = SVector(
                        left_node_back + iy + 1, 
                        left_node_front + iy + 1, 
                        left_node_front + iy
                        )
                    #right faces
                    faces[right_face + (2*iy - 1)] = SVector(
                        left_node_back + m + 2*n - 1 + iy, 
                        left_node_front + m + 2*n - 1 + iy, 
                        left_node_back + m + 2*n - 1 + iy + 1
                        )
                    faces[right_face + (2*iy)] = SVector(
                        left_node_back + m + 2*n - 1 + iy + 1, 
                        left_node_front + m + 2*n - 1 + iy, 
                        left_node_front + m + 2*n - 1 + iy + 1
                        )
                end           
                #ix in n -> right nodes
                nodes[left_node_front + (m + 2*n - 1) + iy] = SVector(
                    a, 
                    (iy - 1)*h, 
                    (iz - 1)*h
                    )
            end
            #iz = p + 1
            if iy != (m + 1)
                #left faces
                faces[4*m*n + 4*n*p + 2*m*p + (p - 1)*2*m + (2*iy - 1)] = SVector(
                    (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + iy, 
                    (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + iy + 1, 
                    (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + iy
                    )
                faces[4*m*n + 4*n*p + 2*m*p + (p - 1)*2*m + (2*iy)] = SVector(
                    (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + iy + 1, 
                    (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + iy + 1, 
                    (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + iy
                    )
                #right faces
                faces[2*m*n + 2*n*p + (p - 1)*2*m + (2*iy - 1)] = SVector(
                    (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 2*n - 1 + iy, 
                    (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + n*(m + 1) + iy, 
                    (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 2*n - 1 + iy + 1
                    )
                faces[2*m*n + 2*n*p + (p - 1)*2*m + (2*iy)] = SVector(
                    (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 2*n - 1 + iy + 1, 
                    (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + n*(m + 1) + iy, 
                    (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + n*(m + 1) + iy + 1
                    )
            end
        end

        # along x - z
        #faces adjacent to the left, right, front and back{
        #iz = 1 {
            #bottom faces
            #ix = 1
            faces[2*m*n + 1] = SVector(1, 
            (n + 1)*(m + 1) + 1, 
            (m + 1) + 1
            )
            faces[2*m*n + 2] = SVector(
                (n + 1)*(m + 1) + 1, 
                (n + 1)*(m + 1) + (m + 1) + 1, 
                (m + 1) + 1
                )
            #ix = n
            faces[2*m*n + (n - 1)*2*p + 1] = SVector(
                (n - 1)*(m + 1) + 1, 
                (n + 1)*(m + 1) + m + 2*n - 2, 
                (n)*(m + 1) + 1
                )
            faces[2*m*n + (n - 1)*2*p + 2] = SVector(
                (n + 1)*(m + 1) + m + 2*n - 2, 
                (n + 1)*(m + 1) + m + 2*n, 
                (n)*(m + 1) + 1
                )
            #top faces
            #ix = 1
            faces[4*m*n + 2*m*p + 2*n*p + 1] = SVector(
                m + 1, 
                2*(m + 1), 
                (n + 1)*(m + 1) + m + 1
                )
            faces[4*m*n + 2*m*p + 2*n*p + 2] = SVector(
                (n + 1)*(m + 1) + m + 1, 
                2*(m + 1), 
                (n + 1)*(m + 1) + m + 3
                )
            #ix = n
            faces[4*m*n + 2*m*p + 2*n*p + (n - 1)*2*p + 1] = SVector(
                n*(m + 1), 
                (n + 1)*(m + 1), 
                (n + 1)*(m + 1) + (m + 1) + 2*(n - 1)
                )
            faces[4*m*n + 2*m*p + 2*n*p + (n - 1)*2*p + 2] = SVector(
                (n + 1)*(m + 1) + (m + 1) + 2*(n - 1), 
                (n + 1)*(m + 1), 
                (n + 1)*(m + 1) + 2*(m + n)
                )
            # }
            #iz = p {
            #ix = 1
            #bottom faces
            faces[2*m*n + (2*p - 1)] = SVector(
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + 1, 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + 1, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 2
                )
            faces[2*m*n + (2*p)] = SVector(
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + 1, 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + m + 2, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 2
                )
            #top faces
            faces[4*m*n + 2*m*p + 2*n*p + (2*p - 1)] = SVector(
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 1, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 3, 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + m + 1
                )
            faces[4*m*n + 2*m*p + 2*n*p + (2*p)] = SVector(
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + m + 1, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + m + 3, 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + 2*(m + 1)
                )
            #ix = n
            #bottom faces
            faces[2*m*n + (n - 1)*2*p + (2*p - 1)] = SVector(
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + (m + 2*(n - 1)), 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + (n - 1)*(m + 1) + 1, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + (m + 2*n)
                )
            faces[2*m*n + (n - 1)*2*p + (2*p)] = SVector(
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + (n - 1)*(m + 1) + 1, 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + (n)*(m + 1) + 1, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + (m + 2*n)
                )
            #top faces
            faces[4*m*n + 2*m*p + 2*n*p + (n - 1)*2*p + (2*p - 1)] = SVector(
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + (m + 2*(n - 1)) + 1, 
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + 2*(m + n), 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + (n)*(m + 1)
                )
            faces[4*m*n + 2*m*p + 2*n*p + (n - 1)*2*p + (2*p)] = SVector(
                (n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + 2*(m + n), 
                2*(n + 1)*(m + 1) + (p - 1)*(2*(m + n)), 
                (n + 1)*(m + 1) + (p - 1)*(2*(m + n)) + (n)*(m + 1)
                )
            # }
            #}
        for ix in range(1, n)
            if (ix != 1) 
                #nodes for iz = p
                #iy = 1 -> bottom nodes
                nodes[(n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + (m + 1) + (ix - 2)*2 + 1] = SVector(
                    (ix - 1)*h, F(0), (p - 1)*h)
                #iy = m -> top nodes
                nodes[(n + 1)*(m + 1) + (p - 2)*(2*(m + n)) + (m + 1) + (ix - 1)*2] = SVector(
                    (ix - 1)*h, b, (p - 1)*h)
                if (ix != n)
                    # iz = 1
                    #bottom faces
                    faces[2*m*n + (ix - 1)*2*p + 1] = SVector(
                        (ix - 1)*(m + 1) + 1, 
                        (n + 1)*(m + 1) + (m + 1) + (ix - 2)*2 + 1, 
                        (ix)*(m + 1) + 1
                        )
                    faces[2*m*n + (ix - 1)*2*p + 2] = SVector(
                        (n + 1)*(m + 1) + (m + 1) + (ix - 2)*2 + 1, 
                        (n + 1)*(m + 1) + (m + 1) + (ix - 1)*2 + 1, 
                        (ix)*(m + 1) + 1
                        )
                    #top faces
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + 1] = SVector(
                        ix*(m + 1), 
                        (ix + 1)*(m + 1), 
                        (n + 1)*(m + 1) + (m + 1) + (ix - 1)*2
                        )
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + 2] = SVector(
                        (n + 1)*(m + 1) + (m + 1) + (ix - 1)*2, 
                        (ix + 1)*(m + 1), 
                        (n + 1)*(m + 1) + (m + 1) + (ix)*2
                        )
                    # iz = p 
                    #bottom faces
                    faces[2*m*n + (ix - 1)*2*p + (2*p - 1)] = SVector(
                        (n + 1)*(m + 1) + (p - 2)*2*(m + n) + (m + 1) + (ix - 2)*2 + 1, 
                        (n + 1)*(m + 1) + (p - 1)*2*(m + n) + (ix - 1)*(m + 1) + 1, 
                        (n + 1)*(m + 1) + (p - 2)*2*(m + n) + (m + 1) + (ix - 1)*2 + 1
                        )
                    faces[2*m*n + (ix - 1)*2*p + (2*p)] = SVector(
                        (n + 1)*(m + 1) + (p - 1)*2*(m + n) + (ix - 1)*(m + 1) + 1, 
                        (n + 1)*(m + 1) + (p - 1)*2*(m + n) + (ix)*(m + 1) + 1, 
                        (n + 1)*(m + 1) + (p - 2)*2*(m + n) + (m + 1) + (ix - 1)*2 + 1
                        )
                    #top faces
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*p - 1)] = SVector(
                        (n + 1)*(m + 1) + (p - 2)*2*(m + n) + (m + 1) + (ix - 1)*2, 
                        (n + 1)*(m + 1) + (p - 2)*2*(m + n) + (m + 1) + (ix)*2, 
                        (n + 1)*(m + 1) + (p - 1)*2*(m + n) + ix*(m + 1)
                        )
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*p)] = SVector(
                        (n + 1)*(m + 1) + (p - 2)*2*(m + n) + (m + 1) + (ix)*2, 
                        (n + 1)*(m + 1) + (p - 1)*2*(m + n) + (ix + 1)*(m + 1), 
                        (n + 1)*(m + 1) + (p - 1)*2*(m + n) + ix*(m + 1)
                        )
                end
            end
            for iz in range(2, p - 1)
                if ix == 1
                    #bottom faces
                    faces[2*m*n + (ix - 1)*2*p + (2*iz - 1)] = SVector(
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + 1, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + 1, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 2
                        )
                    faces[2*m*n + (ix - 1)*2*p + (2*iz)] = SVector(
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + 1, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 2, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 2
                        )
                    #top faces
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*iz - 1)] = SVector(
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 1, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 3, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 1
                        )
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*iz)] = SVector(
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 1, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 3, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 3
                        )
                elseif ix == n
                    #bottom faces
                    faces[2*m*n + (ix - 1)*2*p + (2*iz - 1)] = SVector(
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 2*(n - 1), 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 2*(n - 1), 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 2*n
                        )
                    faces[2*m*n + (ix - 1)*2*p + (2*iz)] = SVector(
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 2*(n - 1), 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 2*n, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 2*n
                        )
                    #top faces
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*iz - 1)] = SVector(
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + m + 1 + 2*(n - 1), 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)), 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 2*(n - 1) + 1
                        )
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*iz)] = SVector(
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)), 
                        (n + 1)*(m + 1) + (iz)*(2*(m + n)), 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + m + 2*(n - 1) + 1
                        )
                    #nodes
                    #iy = 1 -> bottom nodes
                    nodes[(n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 2)*2 + 1] = SVector(
                        (ix - 1)*h, F(0), (iz - 1)*h)
                    #iy = m -> top nodes
                    nodes[(n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 1)*2] = SVector(
                        (ix - 1)*h, b, (iz - 1)*h)
                else
                    #for iy = 1 -> bottom faces
                    faces[2*m*n + (ix - 1)*2*p + (2*iz - 1)] = SVector(
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 2)*(2) + 1, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + (m + 1) + (ix - 2)*(2) + 1, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 1)*2 + 1
                        )
                    faces[2*m*n + (ix - 1)*2*p + (2*iz)] = SVector(
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + (m + 1) + (ix - 2)*(2) + 1, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + (m + 1) + (ix - 1)*(2) + 1, 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 1)*2 + 1
                        )
                    #for iy = m -> top faces
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*iz - 1)] = SVector(
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 1)*(2), 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix)*2, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + (m + 1) + (ix - 1)*(2)
                        )
                    faces[4*m*n + 2*m*p + 2*n*p + (ix - 1)*2*p + (2*iz)] = SVector(
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + (m + 1) + (ix - 1)*(2), 
                        (n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix)*2, 
                        (n + 1)*(m + 1) + (iz - 1)*(2*(m + n)) + (m + 1) + (ix)*(2)
                        )
                    #nodes
                    #iy = 1 -> bottom nodes
                    nodes[(n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 2)*2 + 1] = SVector(
                        (ix - 1)*h, F(0), (iz - 1)*h)
                    #iy = m -> top nodes
                    nodes[(n + 1)*(m + 1) + (iz - 2)*(2*(m + n)) + (m + 1) + (ix - 1)*2] = SVector(
                        (ix - 1)*h, b, (iz - 1)*h)
                end
            end
        end
    return Mesh(nodes, faces)
end

"""
    gmshcuboid(width, height, length, delta)

Creates a mesh for a cuboid of width (along the x-axis) `width`, height (along
    the y-axis) `height` and length (along the z-axis) `length` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

The target edge size is `delta`.
physical => in {"TopPlate", "BottomPlate", "SidePlates", "OpenBox"} extracts and
returns only the specified part of the cuboid

"""
function gmshcuboid(width::T, height::T, length::T, delta::T;physical="ClosedBox", tempname=tempname()) where T
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

//classify parts of geometry
Physical Surface("TopPlate") = {3};
Physical Surface("BottomPlate") = {5};
Physical Surface("SidePlates") = {1,2,4,6};
Physical Surface("OpenBox") = {1,2,4,5,6};
Physical Surface("ClosedBox") = {1,2,3,4,5,6};

Surface Loop(1)={1,2,3,4,5,6};
Volume(1)={1};

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

    # m = read_gmsh_mesh(fno)
    m = read_gmsh_mesh(fno,physical=physical,T=T)

    rm(fno)
    rm(fn)

    return m

end
