# This file is created by reverse engineering Gmsh msh file
# and using ChatGPT for analyzing the output.
# Up to order 5 all results are manually verified.

const gmsh_line_index_to_tuple_order1 = [
    # Vertices
    SVector(1, 0), # V1
    SVector(0, 1), # V2
]

const gmsh_line_index_to_tuple_order2 = [
    # Vertices
    SVector(2, 0), # V1
    SVector(0, 2), # V2

    # Edge
    SVector(1, 1),
]

const gmsh_line_index_to_tuple_order3 = [
    # Vertices
    SVector(3, 0), # V1
    SVector(0, 3), # V2

    # Edge
    SVector(2, 1),
    SVector(1, 2),
]

const gmsh_line_index_to_tuple_order4 = [
    # Vertices
    SVector(4, 0), # V1
    SVector(0, 4), # V2

    # Edge
    SVector(3, 1),
    SVector(2, 2),
    SVector(1, 3),
]

const gmsh_line_index_to_tuple_order5 = [
    # Vertices
    SVector(5, 0), # V1
    SVector(0, 5), # V2

    # Edge
    SVector(4, 1),
    SVector(3, 2),
    SVector(2, 3),
    SVector(1, 4),
]

const gmsh_line_index_to_tuple_order6 = [
    # Vertices
    SVector(6, 0), # V1
    SVector(0, 6), # V2

    # Edge
    SVector(5, 1),
    SVector(4, 2),
    SVector(3, 3),
    SVector(2, 4),
    SVector(1, 5),
]

const gmsh_line_index_to_tuple_order7 = [
    # Vertices
    SVector(7, 0), # V1
    SVector(0, 7), # V2

    # Edge
    SVector(6, 1),
    SVector(5, 2),
    SVector(4, 3),
    SVector(3, 4),
    SVector(2, 5),
    SVector(1, 6),
]

const gmsh_line_index_to_tuple_order8 = [
    # Vertices
    SVector(8, 0), # V1
    SVector(0, 8), # V2

    # Edge
    SVector(7, 1),
    SVector(6, 2),
    SVector(5, 3),
    SVector(4, 4),
    SVector(3, 5),
    SVector(2, 6),
    SVector(1, 7),
]

const gmsh_line_index_to_tuple_order9 = [
    # Vertices
    SVector(9, 0), # V1
    SVector(0, 9), # V2

    # Edge
    SVector(8, 1),
    SVector(7, 2),
    SVector(6, 3),
    SVector(5, 4),
    SVector(4, 5),
    SVector(3, 6),
    SVector(2, 7),
    SVector(1, 8),
]

const gmsh_line_index_to_tuple_order10 = [
    # Vertices
    SVector(10, 0), # V1
    SVector(0, 10), # V2

    # Edge
    SVector(9, 1),
    SVector(8, 2),
    SVector(7, 3),
    SVector(6, 4),
    SVector(5, 5),
    SVector(4, 6),
    SVector(3, 7),
    SVector(2, 8),
    SVector(1, 9),
]

const gmsh_line_index_to_tuple = [
    gmsh_line_index_to_tuple_order1,
    gmsh_line_index_to_tuple_order2,
    gmsh_line_index_to_tuple_order3,
    gmsh_line_index_to_tuple_order4,
    gmsh_line_index_to_tuple_order5,
    gmsh_line_index_to_tuple_order6,
    gmsh_line_index_to_tuple_order7,
    gmsh_line_index_to_tuple_order8,
    gmsh_line_index_to_tuple_order9,
    gmsh_line_index_to_tuple_order10,
]