using CompScienceMeshes
using Gmsh
using StaticArrays
using Test

# Get all nodes
gmsh.initialize()
gmsh.open(joinpath(@__DIR__, "mesa.msh"))

node_tags, node_coords, _ = gmsh.model.mesh.getNodes(
    -1, # dim<0: Get nodes of entities of any dimension (all node)
    -1, # tag<0: Get the nods of all entities of dimension dim
    true, # includeBoundary (default false)
    true #  returnParametricCoord (default true)
)

element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements()

gmsh.finalize()

vertices3d =  CompScienceMeshes.gmsh_to_mesh_nodes(node_tags, node_coords)
vertices2d =  CompScienceMeshes.gmsh_to_mesh_nodes(node_tags, node_coords, udim=2)
vertices3df =  CompScienceMeshes.gmsh_to_mesh_nodes(node_tags, node_coords, vertextype=Float32)
vertices2df =  CompScienceMeshes.gmsh_to_mesh_nodes(
    node_tags,
    node_coords,
    udim=2,
    vertextype=Float32
)

@test eltype(first(vertices3d)) == Float64
@test length(first(vertices3d)) == 3
@test eltype(first(vertices2d)) == Float64
@test length(first(vertices2d)) == 2
@test eltype(first(vertices3df)) == Float32
@test length(first(vertices3df)) == 3
@test eltype(first(vertices2df)) == Float32
@test length(first(vertices2df)) == 2

elements = CompScienceMeshes.gmsh_get_elements(element_types,
    element_tags,
    element_node_tags;
    element=:triangle,
    order=1
)
@test length(elements) == length(element_tags[2])

# Ensure backwards compatibility: Triangular surface

mesh_old_interface = CompScienceMeshes.read_gmsh_mesh(joinpath(@__DIR__, "mesa.msh"))
mesh_new_interface = load_gmsh_mesh(joinpath(@__DIR__, "mesa.msh"))

@test mesh_new_interface.vertices == mesh_old_interface.vertices
@test mesh_new_interface.faces == mesh_old_interface.faces

# Ensure backwards compatibility: Trianglar surface with physical subdomain

mesh_old_interface = CompScienceMeshes.read_gmsh_mesh(
    joinpath(@__DIR__, "assets/cuboid_triangle.msh"),
    physical="Top"
)
mesh_new_interface = load_gmsh_mesh(
    joinpath(@__DIR__, "assets/cuboid_triangle.msh"),
    physical="Top"
)

@test mesh_new_interface.vertices == mesh_old_interface.vertices
@test mesh_new_interface.faces == mesh_old_interface.faces

# Ensure backwards compatibility: Tetrahedra

mesh_old_interface = CompScienceMeshes.read_gmsh3d_mesh(
    joinpath(@__DIR__, "assets/cuboid_tetrahedron.msh")
)
mesh_new_interface = load_gmsh_mesh(
    joinpath(@__DIR__, "assets/cuboid_tetrahedron.msh"),
    element=:tetrahedron,
    sort=false
)

@test mesh_new_interface.vertices == mesh_old_interface.vertices

# TODO: read_gmsh3d_mesh performs a reordering of the vertices of the tetrahedrons
# to ensure that they are all oriented likewise
# I am not sure this is still needed.
# In particular, the higher-order case would require a careful treatment

# @test mesh_new_interface.faces == mesh_old_interface.faces 

@test CompScienceMeshes.isoriented(mesh_new_interface) == true

# Test loading of quadrangle mesh

mesh_new_interface = load_gmsh_mesh(
    joinpath(@__DIR__, "assets/cuboid_quadrangle.msh"),
    element=:quadrangle
)

# Ensure backwards compatibility: Triangles + Physical Domain

mesh_old_interface = CompScienceMeshes.read_gmsh_mesh(
    joinpath(@__DIR__, "assets/cuboid_triangle.msh"),
    physical="Top"
)
mesh_new_interface = load_gmsh_mesh(
    joinpath(@__DIR__, "assets/cuboid_triangle.msh"),
    physical="Top"
)

@test mesh_new_interface.vertices == mesh_old_interface.vertices
@test mesh_new_interface.faces == mesh_old_interface.faces

# Ensure backwards compatibility: Tetrahedrons + Physical Domain

mesh_old_interface = CompScienceMeshes.read_gmsh3d_mesh(
    joinpath(@__DIR__, "assets/cuboid_tetrahedron.msh"),
    physical="Interior"
)
mesh_new_interface = load_gmsh_mesh(joinpath(@__DIR__, "assets/cuboid_tetrahedron.msh"),
    element=:tetrahedron,
    physical="Interior"
)

@test mesh_new_interface.vertices == mesh_old_interface.vertices
#@test mesh_new_interface.faces == mesh_old_interface.faces
@test CompScienceMeshes.isoriented(mesh_new_interface) == true

# Load 2D geometry
mesh2d = load_gmsh_mesh(
    joinpath(@__DIR__, "assets/circle2d.msh"),
    udim=2,
    vertextype=Float32,
    element=:line
)

@test typeof(vertices(mesh2d)) == Vector{SVector{2, Float32}}

# Load 3D geometry with lines
mesh3d = load_gmsh_mesh(
    joinpath(@__DIR__, "assets/circle2d.msh"),
    udim=3,
    vertextype=Float32,
    element=:line
)

@test typeof(vertices(mesh2d)) == Vector{SVector{2, Float32}}