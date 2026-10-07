# This file is created by reverse engineering Gmsh msh file
# and using ChatGPT for analyzing the output.
# Up to order 5 all results are manually verified.

const gmsh_triangle_index_to_tuple_order1 = [
    # Vertices
    SVector(1, 0, 0), # V1
    SVector(0, 1, 0), # V2
    SVector(0, 0, 1), # V3
]

const gmsh_triangle_index_to_tuple_order2 = [
    # Vertices
    SVector(2, 0, 0), # V1
    SVector(0, 2, 0), # V2
    SVector(0, 0, 2), # V3

    # Edge 1: V1 -> V2
    SVector(1, 1, 0),

    # Edge 2: V2 -> V3
    SVector(0, 1, 1),

    # Edge 3: V2 -> V3
    SVector(1, 0, 1),
]

const gmsh_triangle_index_to_tuple_order3 = [
    # Vertices
    SVector(3, 0, 0), # V1
    SVector(0, 3, 0), # V2
    SVector(0, 0, 3), # V3

    # Edge 1: V1 -> V2
    SVector(2, 1, 0),
    SVector(1, 2, 0),

    # Edge 2: V2 -> V3
    SVector(0, 2, 1),
    SVector(0, 1, 2),

    # Edge 3: V2 -> V3
    SVector(1, 0, 2),
    SVector(2, 0, 1),

    # Interior
    SVector(1, 1, 1),
]

const gmsh_triangle_index_to_tuple_order4 = [
    # Vertices
    SVector(4, 0, 0), # V1
    SVector(0, 4, 0), # V2
    SVector(0, 0, 4), # V3

    # Edge 1: V1 -> V2
    SVector(3, 1, 0),
    SVector(2, 2, 0),
    SVector(1, 3, 0),

    # Edge 2: V2 -> V3
    SVector(0, 3, 1),
    SVector(0, 2, 2),
    SVector(0, 1, 3),

    # Edge 3: V3 -> V1
    SVector(1, 0, 3),
    SVector(2, 0, 2),
    SVector(3, 0, 1),

    # Interior
    SVector(2, 1, 1),
    SVector(1, 2, 1),
    SVector(1, 1, 2),
]

const gmsh_triangle_index_to_tuple_order5 = [
    # Vertices
    SVector(5, 0, 0), # V1
    SVector(0, 5, 0), # V2
    SVector(0, 0, 5), # V3

    # Edge 1: V1 -> V2
    SVector(4, 1, 0),
    SVector(3, 2, 0),
    SVector(2, 3, 0),
    SVector(1, 4, 0),

    # Edge 2: V2 -> V3
    SVector(0, 4, 1),
    SVector(0, 3, 2),
    SVector(0, 2, 3),
    SVector(0, 1, 4),

    # Edge 3: V3 -> V1
    SVector(1, 0, 4),
    SVector(2, 0, 3),
    SVector(3, 0, 2),
    SVector(4, 0, 1),

    # Interior
    SVector(3, 1, 1),
    SVector(1, 3, 1),
    SVector(1, 1, 3),
    SVector(2, 2, 1),
    SVector(1, 2, 2),
    SVector(2, 1, 2),
]

const gmsh_triangle_index_to_tuple_order6 = [
    # Vertices
    SVector(6, 0, 0), # V1
    SVector(0, 6, 0), # V2
    SVector(0, 0, 6), # V3

    # Edge 1: V1 -> V2
    SVector(5, 1, 0),
    SVector(4, 2, 0),
    SVector(3, 3, 0),
    SVector(2, 4, 0),
    SVector(1, 5, 0),

    # Edge 2: V2 -> V3
    SVector(0, 5, 1),
    SVector(0, 4, 2),
    SVector(0, 3, 3),
    SVector(0, 2, 4),
    SVector(0, 1, 5),

    # Edge 3: V3 -> V1
    SVector(1, 0, 5),
    SVector(2, 0, 4),
    SVector(3, 0, 3),
    SVector(4, 0, 2),
    SVector(5, 0, 1),

    # Interior nodes
    SVector(4, 1, 1),
    SVector(1, 4, 1),
    SVector(1, 1, 4),
    SVector(3, 2, 1),
    SVector(2, 3, 1),
    SVector(1, 3, 2),
    SVector(1, 2, 3),
    SVector(2, 1, 3),
    SVector(3, 1, 2),
    SVector(2, 2, 2),
]

const gmsh_triangle_index_to_tuple_order7 = [
    # Vertices
    SVector(7, 0, 0), # V1
    SVector(0, 7, 0), # V2
    SVector(0, 0, 7), # V3

    # Edge 1: V1 -> V2
    SVector(6, 1, 0),
    SVector(5, 2, 0),
    SVector(4, 3, 0),
    SVector(3, 4, 0),
    SVector(2, 5, 0),
    SVector(1, 6, 0),

    # Edge 2: V2 -> V3
    SVector(0, 6, 1),
    SVector(0, 5, 2),
    SVector(0, 4, 3),
    SVector(0, 3, 4),
    SVector(0, 2, 5),
    SVector(0, 1, 6),

    # Edge 3: V3 -> V1
    SVector(1, 0, 6),
    SVector(2, 0, 5),
    SVector(3, 0, 4),
    SVector(4, 0, 3),
    SVector(5, 0, 2),
    SVector(6, 0, 1),

    # Interior
    SVector(5, 1, 1),
    SVector(1, 5, 1),
    SVector(1, 1, 5),
    SVector(4, 2, 1),
    SVector(3, 3, 1),
    SVector(2, 4, 1),
    SVector(1, 4, 2),
    SVector(1, 3, 3),
    SVector(1, 2, 4),
    SVector(2, 1, 4),
    SVector(3, 1, 3),
    SVector(4, 1, 2),
    SVector(3, 2, 2),
    SVector(2, 3, 2),
    SVector(2, 2, 3),
]

const gmsh_triangle_index_to_tuple_order8 = [
    # Vertices
    SVector(8, 0, 0), # V1
    SVector(0, 8, 0), # V2
    SVector(0, 0, 8), # V3

    # Edge 1: V1 -> V2
    SVector(7, 1, 0),
    SVector(6, 2, 0),
    SVector(5, 3, 0),
    SVector(4, 4, 0),
    SVector(3, 5, 0),
    SVector(2, 6, 0),
    SVector(1, 7, 0),

    # Edge 2: V2 -> V3
    SVector(0, 7, 1),
    SVector(0, 6, 2),
    SVector(0, 5, 3),
    SVector(0, 4, 4),
    SVector(0, 3, 5),
    SVector(0, 2, 6),
    SVector(0, 1, 7),

    # Edge 3: V3 -> V1
    SVector(1, 0, 7),
    SVector(2, 0, 6),
    SVector(3, 0, 5),
    SVector(4, 0, 4),
    SVector(5, 0, 3),
    SVector(6, 0, 2),
    SVector(7, 0, 1),

    # Interior
    SVector(6, 1, 1),
    SVector(1, 6, 1),
    SVector(1, 1, 6),
    SVector(5, 2, 1),
    SVector(4, 3, 1),
    SVector(3, 4, 1),
    SVector(2, 5, 1),
    SVector(1, 5, 2),
    SVector(1, 4, 3),
    SVector(1, 3, 4),
    SVector(1, 2, 5),
    SVector(2, 1, 5),
    SVector(3, 1, 4),
    SVector(4, 1, 3),
    SVector(5, 1, 2),
    SVector(4, 2, 2),
    SVector(2, 4, 2),
    SVector(2, 2, 4),
    SVector(3, 3, 2),
    SVector(2, 3, 3),
    SVector(3, 2, 3),
]

const gmsh_triangle_index_to_tuple_order9 = [
    # Vertices
    SVector(9, 0, 0), # V1
    SVector(0, 9, 0), # V2
    SVector(0, 0, 9), # V3

    # Edge 1: V1 -> V2
    SVector(8, 1, 0),
    SVector(7, 2, 0),
    SVector(6, 3, 0),
    SVector(5, 4, 0),
    SVector(4, 5, 0),
    SVector(3, 6, 0),
    SVector(2, 7, 0),
    SVector(1, 8, 0),

    # Edge 2: V2 -> V3
    SVector(0, 8, 1),
    SVector(0, 7, 2),
    SVector(0, 6, 3),
    SVector(0, 5, 4),
    SVector(0, 4, 5),
    SVector(0, 3, 6),
    SVector(0, 2, 7),
    SVector(0, 1, 8),

    # Edge 3: V3 -> V1
    SVector(1, 0, 8),
    SVector(2, 0, 7),
    SVector(3, 0, 6),
    SVector(4, 0, 5),
    SVector(5, 0, 4),
    SVector(6, 0, 3),
    SVector(7, 0, 2),
    SVector(8, 0, 1),

    # Interior
    SVector(7, 1, 1),
    SVector(1, 7, 1),
    SVector(1, 1, 7),
    SVector(6, 2, 1),
    SVector(5, 3, 1),
    SVector(4, 4, 1),
    SVector(3, 5, 1),
    SVector(2, 6, 1),
    SVector(1, 6, 2),
    SVector(1, 5, 3),
    SVector(1, 4, 4),
    SVector(1, 3, 5),
    SVector(1, 2, 6),
    SVector(2, 1, 6),
    SVector(3, 1, 5),
    SVector(4, 1, 4),
    SVector(5, 1, 3),
    SVector(6, 1, 2),
    SVector(5, 2, 2),
    SVector(2, 5, 2),
    SVector(2, 2, 5),
    SVector(4, 3, 2),
    SVector(3, 4, 2),
    SVector(2, 4, 3),
    SVector(2, 3, 4),
    SVector(3, 2, 4),
    SVector(4, 2, 3),
    SVector(3, 3, 3),
]

const gmsh_triangle_index_to_tuple_order10 = [
    # Vertices
    SVector(10, 0, 0),
    SVector(0, 10, 0),
    SVector(0, 0, 10),

    # Edge 1
    SVector(9, 1, 0),
    SVector(8, 2, 0),
    SVector(7, 3, 0),
    SVector(6, 4, 0),
    SVector(5, 5, 0),
    SVector(4, 6, 0),
    SVector(3, 7, 0),
    SVector(2, 8, 0),
    SVector(1, 9, 0),

    # Edge 2
    SVector(0, 9, 1),
    SVector(0, 8, 2),
    SVector(0, 7, 3),
    SVector(0, 6, 4),
    SVector(0, 5, 5),
    SVector(0, 4, 6),
    SVector(0, 3, 7),
    SVector(0, 2, 8),
    SVector(0, 1, 9),

    # Edge 3
    SVector(1, 0, 9),
    SVector(2, 0, 8),
    SVector(3, 0, 7),
    SVector(4, 0, 6),
    SVector(5, 0, 5),
    SVector(6, 0, 4),
    SVector(7, 0, 3),
    SVector(8, 0, 2),
    SVector(9, 0, 1),

    # Interior nodes
    SVector(8, 1, 1),
    SVector(1, 8, 1),
    SVector(1, 1, 8),
    SVector(7, 2, 1),
    SVector(6, 3, 1),
    SVector(5, 4, 1),
    SVector(4, 5, 1),
    SVector(3, 6, 1),
    SVector(2, 7, 1),
    SVector(1, 7, 2),
    SVector(1, 6, 3),
    SVector(1, 5, 4),
    SVector(1, 4, 5),
    SVector(1, 3, 6),
    SVector(1, 2, 7),
    SVector(2, 1, 7),
    SVector(3, 1, 6),
    SVector(4, 1, 5),
    SVector(5, 1, 4),
    SVector(6, 1, 3),
    SVector(7, 1, 2),
    SVector(6, 2, 2),
    SVector(2, 6, 2),
    SVector(2, 2, 6),
    SVector(5, 3, 2),
    SVector(4, 4, 2),
    SVector(3, 5, 2),
    SVector(2, 5, 3),
    SVector(2, 4, 4),
    SVector(2, 3, 5),
    SVector(3, 2, 5),
    SVector(4, 2, 4),
    SVector(5, 2, 3),
    SVector(4, 3, 3),
    SVector(3, 4, 3),
    SVector(3, 3, 4),
]

const gmsh_triangle_index_to_tuple = [
    gmsh_triangle_index_to_tuple_order1,
    gmsh_triangle_index_to_tuple_order2,
    gmsh_triangle_index_to_tuple_order3,
    gmsh_triangle_index_to_tuple_order4,
    gmsh_triangle_index_to_tuple_order5,
    gmsh_triangle_index_to_tuple_order6,
    gmsh_triangle_index_to_tuple_order7,
    gmsh_triangle_index_to_tuple_order8,
    gmsh_triangle_index_to_tuple_order9,
    gmsh_triangle_index_to_tuple_order10,
]