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
        for ix in 1:n
            for iy in 1:m
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
        for iy in 0:m
            #front
            nodes[n*(m + 1) + iy + 1] = SVector(n*h, (iy*h),  F(0))
            #back
            nodes[back_node + n*(m + 1) + iy + 1] = SVector(
                n*h, (iy*h), p*h) 
        end 

        # along y - z
        for iy in 1 : m + 1
            for iz in 2 : p
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
        for ix in 1 : n
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
            for iz in 2 : p - 1
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



"""
    meshcuboid1hole(width, height, holewidth, h)

Create a mesh of a cuboid of size `width × width × height` with a hole of size `holewidth × holewidth × height` (̂x × ̂y × ̂z, resp.).

The target edge size is `h`.
"""
function meshcuboid1hole(width, height, holewidth, h)
    @assert holewidth < width
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("squaretorus")

    # bottom plate
    gmsh.model.geo.addPoint(holewidth/2, -holewidth/2, -height/2, h, 1)
    gmsh.model.geo.addPoint(holewidth/2, holewidth/2, -height/2, h, 2)
    gmsh.model.geo.addPoint(-holewidth/2, holewidth/2, -height/2, h, 3)
    gmsh.model.geo.addPoint(-holewidth/2, -holewidth/2, -height/2, h, 4)
    gmsh.model.geo.addPoint(width/2, -width/2, -height/2, h, 5)
    gmsh.model.geo.addPoint(width/2, width/2, -height/2, h, 6)
    gmsh.model.geo.addPoint(-width/2, width/2, -height/2, h, 7)
    gmsh.model.geo.addPoint(-width/2, -width/2, -height/2, h, 8)

    gmsh.model.geo.addLine(2, 3, 1)
    gmsh.model.geo.addLine(3, 4, 2)
    gmsh.model.geo.addLine(4, 1, 3)
    gmsh.model.geo.addLine(1, 2, 4)
    gmsh.model.geo.addCurveLoop([-1, -2, -3, -4], 101)

    gmsh.model.geo.addLine(6, 7, 5)
    gmsh.model.geo.addLine(7, 8, 6)
    gmsh.model.geo.addLine(8, 5, 7)
    gmsh.model.geo.addLine(5, 6, 8)
    gmsh.model.geo.addCurveLoop([-5, -6, -7, -8], 102)
    gmsh.model.geo.addPlaneSurface([-101, -102], 1)

    # top plate
    gmsh.model.geo.addPoint(holewidth/2, -holewidth/2, height/2, h, 9)
    gmsh.model.geo.addPoint(holewidth/2, holewidth/2, height/2, h, 10)
    gmsh.model.geo.addPoint(-holewidth/2, holewidth/2, height/2, h, 11)
    gmsh.model.geo.addPoint(-holewidth/2, -holewidth/2, height/2, h, 12)
    gmsh.model.geo.addPoint(width/2, -width/2, height/2, h, 13)
    gmsh.model.geo.addPoint(width/2, width/2, height/2, h, 14)
    gmsh.model.geo.addPoint(-width/2, width/2, height/2, h, 15)
    gmsh.model.geo.addPoint(-width/2, -width/2, height/2, h, 16)

    gmsh.model.geo.addLine(10, 11, 9)
    gmsh.model.geo.addLine(11, 12, 10)
    gmsh.model.geo.addLine(12, 9, 11)
    gmsh.model.geo.addLine(9, 10, 12)
    gmsh.model.geo.addCurveLoop([9, 10, 11, 12], 103)
    
    gmsh.model.geo.addLine(14, 15, 13)
    gmsh.model.geo.addLine(15, 16, 14)
    gmsh.model.geo.addLine(16, 13, 15)
    gmsh.model.geo.addLine(13, 14, 16)
    gmsh.model.geo.addCurveLoop([13, 14, 15, 16], 104)
    gmsh.model.geo.addPlaneSurface([-103, -104], 2)

    # sides
    gmsh.model.geo.addLine(2, 10, 17)
    gmsh.model.geo.addLine(3, 11, 18)
    gmsh.model.geo.addLine(4, 12, 19)
    gmsh.model.geo.addLine(1, 9, 20)
    gmsh.model.geo.addLine(6, 14, 21)
    gmsh.model.geo.addLine(7, 15, 22)
    gmsh.model.geo.addLine(8, 16, 23)
    gmsh.model.geo.addLine(5, 13, 24)

    gmsh.model.geo.addCurveLoop([-12, 17, 4, -20], 105)
    gmsh.model.geo.addPlaneSurface([-105], 3)

    gmsh.model.geo.addCurveLoop([1, 18, -9, -17], 106)
    gmsh.model.geo.addPlaneSurface([-106], 4)

    gmsh.model.geo.addCurveLoop([-18, -10, 19, 2], 107)
    gmsh.model.geo.addPlaneSurface([-107], 5)

    gmsh.model.geo.addCurveLoop([-19, -11, 20, 3], 108)
    gmsh.model.geo.addPlaneSurface([-108], 6)

    gmsh.model.geo.addCurveLoop([24, 16, -21, -8], 109)
    gmsh.model.geo.addPlaneSurface([-109], 7)

    gmsh.model.geo.addCurveLoop([-5, -22, 13, 21], 110)
    gmsh.model.geo.addPlaneSurface([-110], 8)

    gmsh.model.geo.addCurveLoop([-6, -23, 14, 22], 111)
    gmsh.model.geo.addPlaneSurface([-111], 9)

    gmsh.model.geo.addCurveLoop([23, 15, -24, -7], 112)
    gmsh.model.geo.addPlaneSurface([-112], 10)

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
    meshcuboid4holes(width, height, holewidth, h)

Create a mesh of a cuboid of size `width × width × height` with 4 holes of size `holewidth × holewidth × height` (̂x × ̂y × ̂z, resp.).

The target edge size is `h`.
"""

function meshcuboid4holes(width, height, holewidth, h)
    @assert 2*holewidth < width
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("squaretorus4holes")

    # bottom plate
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 - holewidth/2, -height/2, h, 1)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 + holewidth/2, -height/2, h, 2)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 + holewidth/2, -height/2, h, 3)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 - holewidth/2, -height/2, h, 4)

    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 - holewidth/2, -height/2, h, 5)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 + holewidth/2, -height/2, h, 6)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 + holewidth/2, -height/2, h, 7)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 - holewidth/2, -height/2, h, 8)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 - holewidth/2, -height/2, h, 9)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 + holewidth/2, -height/2, h, 10)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 + holewidth/2, -height/2, h, 11)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 - holewidth/2, -height/2, h, 12)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 - holewidth/2, -height/2, h, 13)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 + holewidth/2, -height/2, h, 14)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 + holewidth/2, -height/2, h, 15)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 - holewidth/2, -height/2, h, 16)

    gmsh.model.geo.addPoint(width/2, -width/2, -height/2, h, 17)
    gmsh.model.geo.addPoint(width/2, width/2, -height/2, h, 18)
    gmsh.model.geo.addPoint(-width/2, width/2, -height/2, h, 19)
    gmsh.model.geo.addPoint(-width/2, -width/2, -height/2, h, 20)

    gmsh.model.geo.addLine(2, 3, 1)
    gmsh.model.geo.addLine(3, 4, 2)
    gmsh.model.geo.addLine(4, 1, 3)
    gmsh.model.geo.addLine(1, 2, 4)
    gmsh.model.geo.addCurveLoop([-1, -2, -3, -4], 101)

    gmsh.model.geo.addLine(6, 7, 5)
    gmsh.model.geo.addLine(7, 8, 6)
    gmsh.model.geo.addLine(8, 5, 7)
    gmsh.model.geo.addLine(5, 6, 8)
    gmsh.model.geo.addCurveLoop([-5, -6, -7, -8], 102)

    gmsh.model.geo.addLine(10, 11, 9)
    gmsh.model.geo.addLine(11, 12, 10)
    gmsh.model.geo.addLine(12, 9, 11)
    gmsh.model.geo.addLine(9, 10, 12)
    gmsh.model.geo.addCurveLoop([-9, -10, -11, -12], 103)

    gmsh.model.geo.addLine(14, 15, 13)
    gmsh.model.geo.addLine(15, 16, 14)
    gmsh.model.geo.addLine(16, 13, 15)
    gmsh.model.geo.addLine(13, 14, 16)
    gmsh.model.geo.addCurveLoop([-13, -14, -15, -16], 104)

    gmsh.model.geo.addLine(18, 19, 17)
    gmsh.model.geo.addLine(19, 20, 18)
    gmsh.model.geo.addLine(20, 17, 19)
    gmsh.model.geo.addLine(17, 18, 20)
    gmsh.model.geo.addCurveLoop([-17, -18, -19, -20], 105)

    gmsh.model.geo.addPlaneSurface([-101, -102, -103, -104, -105], 1)

    # top plate
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 - holewidth/2, height/2, h, 21)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, -width/4 + holewidth/2, height/2, h, 22)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 + holewidth/2, height/2, h, 23)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, -width/4 - holewidth/2, height/2, h, 24)

    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 - holewidth/2, height/2, h, 25)
    gmsh.model.geo.addPoint(width/4 + holewidth/2, width/4 + holewidth/2, height/2, h, 26)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 + holewidth/2, height/2, h, 27)
    gmsh.model.geo.addPoint(width/4 - holewidth/2, width/4 - holewidth/2, height/2, h, 28)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 - holewidth/2, height/2, h, 29)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, width/4 + holewidth/2, height/2, h, 30)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 + holewidth/2, height/2, h, 31)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, width/4 - holewidth/2, height/2, h, 32)

    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 - holewidth/2, height/2, h, 33)
    gmsh.model.geo.addPoint(-width/4 + holewidth/2, -width/4 + holewidth/2, height/2, h, 34)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 + holewidth/2, height/2, h, 35)
    gmsh.model.geo.addPoint(-width/4 - holewidth/2, -width/4 - holewidth/2, height/2, h, 36)

    gmsh.model.geo.addPoint(width/2, -width/2, height/2, h, 37)
    gmsh.model.geo.addPoint(width/2, width/2, height/2, h, 38)
    gmsh.model.geo.addPoint(-width/2, width/2, height/2, h, 39)
    gmsh.model.geo.addPoint(-width/2, -width/2, height/2, h, 40)

    gmsh.model.geo.addLine(22, 23, 21)
    gmsh.model.geo.addLine(23, 24, 22)
    gmsh.model.geo.addLine(24, 21, 23)
    gmsh.model.geo.addLine(21, 22, 24)
    gmsh.model.geo.addCurveLoop([-21, -22, -23, -24], 106)

    gmsh.model.geo.addLine(26, 27, 25)
    gmsh.model.geo.addLine(27, 28, 26)
    gmsh.model.geo.addLine(28, 25, 27)
    gmsh.model.geo.addLine(25, 26, 28)
    gmsh.model.geo.addCurveLoop([-25, -26, -27, -28], 107)

    gmsh.model.geo.addLine(30, 31, 29)
    gmsh.model.geo.addLine(31, 32, 30)
    gmsh.model.geo.addLine(32, 29, 31)
    gmsh.model.geo.addLine(29, 30, 32)
    gmsh.model.geo.addCurveLoop([-29, -30, -31, -32], 108)

    gmsh.model.geo.addLine(34, 35, 33)
    gmsh.model.geo.addLine(35, 36, 34)
    gmsh.model.geo.addLine(36, 33, 35)
    gmsh.model.geo.addLine(33, 34, 36)
    gmsh.model.geo.addCurveLoop([-33, -34, -35, -36], 109)

    gmsh.model.geo.addLine(38, 39, 37)
    gmsh.model.geo.addLine(39, 40, 38)
    gmsh.model.geo.addLine(40, 37, 39)
    gmsh.model.geo.addLine(37, 38, 40)
    gmsh.model.geo.addCurveLoop([-37, -38, -39, -40], 110)

    gmsh.model.geo.addPlaneSurface([106, 107, 108, 109, 110], 2)

    # sides
    gmsh.model.geo.addLine(2, 22, 41)
    gmsh.model.geo.addLine(3, 23, 42)
    gmsh.model.geo.addLine(4, 24, 43)
    gmsh.model.geo.addLine(1, 21, 44)
    gmsh.model.geo.addLine(6, 26, 45)
    gmsh.model.geo.addLine(7, 27, 46)
    gmsh.model.geo.addLine(8, 28, 47)
    gmsh.model.geo.addLine(5, 25, 48)
    gmsh.model.geo.addLine(10, 30, 49)
    gmsh.model.geo.addLine(11, 31, 50)
    gmsh.model.geo.addLine(12, 32, 51)
    gmsh.model.geo.addLine(9, 29, 52)
    gmsh.model.geo.addLine(14, 34, 53)
    gmsh.model.geo.addLine(15, 35, 54)
    gmsh.model.geo.addLine(16, 36, 55)
    gmsh.model.geo.addLine(13, 33, 56)
    gmsh.model.geo.addLine(18, 38, 57)
    gmsh.model.geo.addLine(19, 39, 58)
    gmsh.model.geo.addLine(20, 40, 59)
    gmsh.model.geo.addLine(17, 37, 60)

    gmsh.model.geo.addCurveLoop([1, 42, -21, -41], 111)
    gmsh.model.geo.addPlaneSurface([-111], 3)

    gmsh.model.geo.addCurveLoop([2, 43, -22, -42], 112)
    gmsh.model.geo.addPlaneSurface([-112], 4)

    gmsh.model.geo.addCurveLoop([3, 44, -23, -43], 113)
    gmsh.model.geo.addPlaneSurface([-113], 5)

    gmsh.model.geo.addCurveLoop([4, 41, -24, -44], 114)
    gmsh.model.geo.addPlaneSurface([-114], 6)


    gmsh.model.geo.addCurveLoop([5, 46, -25, -45], 115)
    gmsh.model.geo.addPlaneSurface([-115], 7)

    gmsh.model.geo.addCurveLoop([6, 47, -26, -46], 116)
    gmsh.model.geo.addPlaneSurface([-116], 8)

    gmsh.model.geo.addCurveLoop([7, 48, -27, -47], 117)
    gmsh.model.geo.addPlaneSurface([-117], 9)

    gmsh.model.geo.addCurveLoop([8, 45, -28, -48], 118)
    gmsh.model.geo.addPlaneSurface([-118], 10)
    

    gmsh.model.geo.addCurveLoop([9, 50, -29, -49], 119)
    gmsh.model.geo.addPlaneSurface([-119], 11)

    gmsh.model.geo.addCurveLoop([10, 51, -30, -50], 120)
    gmsh.model.geo.addPlaneSurface([-120], 12)

    gmsh.model.geo.addCurveLoop([11, 52, -31, -51], 121)
    gmsh.model.geo.addPlaneSurface([-121], 13)

    gmsh.model.geo.addCurveLoop([12, 49, -32, -52], 122)
    gmsh.model.geo.addPlaneSurface([-122], 14)

    
    gmsh.model.geo.addCurveLoop([13, 54, -33, -53], 123)
    gmsh.model.geo.addPlaneSurface([-123], 15)

    gmsh.model.geo.addCurveLoop([14, 55, -34, -54], 124)
    gmsh.model.geo.addPlaneSurface([-124], 16)

    gmsh.model.geo.addCurveLoop([15, 56, -35, -55], 125)
    gmsh.model.geo.addPlaneSurface([-125], 17)

    gmsh.model.geo.addCurveLoop([16, 53, -36, -56], 126)
    gmsh.model.geo.addPlaneSurface([-126], 18)

    
    gmsh.model.geo.addCurveLoop([17, 58, -37, -57], 127)
    gmsh.model.geo.addPlaneSurface([127], 19)

    gmsh.model.geo.addCurveLoop([18, 59, -38, -58], 128)
    gmsh.model.geo.addPlaneSurface([128], 20)

    gmsh.model.geo.addCurveLoop([19, 60, -39, -59], 129)
    gmsh.model.geo.addPlaneSurface([129], 21)

    gmsh.model.geo.addCurveLoop([20, 57, -40, -60], 130)
    gmsh.model.geo.addPlaneSurface([130], 22)

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