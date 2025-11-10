using Test
using CompScienceMeshes
using StaticArrays

## Line order 1
v1 = SVector(1.0, 0.0)
v2 = SVector(0.0, 0.0)

verts1 = [v1, v2]
faces1 = [SVector(1,2)]

m1 = CompScienceMeshes.Mesh(verts1, faces1)
m2 = CompScienceMeshes.CurvilinearMesh(verts1, faces1, 1, 1)

ch1 = chart(m1, 1)
ch2 = chart(m2, 1)

@test barytocart(ch1, 1.0) == SVector(1.0, 0.0)
@test barytocart(ch2, 1.0) == SVector(1.0, 0.0)
@test barytocart(ch1, 0.0) == SVector(0.0, 0.0)
@test barytocart(ch2, 0.0) == SVector(0.0, 0.0)

mp1 = neighborhood(ch1, 0.0)
mp2 = neighborhood(ch2, 0.0)

@test tangents(mp1, 1) == SVector(1.0, 0.0)
@test tangents(mp2, 1) == SVector(1.0, 0.0)

## Line order 2

v3 = SVector(0.5, 0.0)
verts2 = [v1, v2, v3]
faces2 = [SVector(1,2,3)]

m2 = CompScienceMeshes.CurvilinearMesh(verts2, faces2, 1, 2)

ch2 = chart(m2, 1)

@test barytocart(ch1, 1.0) == SVector(1.0, 0.0)
@test barytocart(ch2, 1.0) == SVector(1.0, 0.0)
@test barytocart(ch1, 0.0) == SVector(0.0, 0.0)
@test barytocart(ch2, 0.0) == SVector(0.0, 0.0)

mp1 = neighborhood(ch1, 0.0)
mp2 = neighborhood(ch2, 0.0)

@test tangents(mp1, 1) == SVector(1.0, 0.0)
@test tangents(mp2, 1) == SVector(1.0, 0.0) # Sign is flipped

## Line order 2  (parabola)

v1 = SVector(0.0, 0.0^2)
v2 = SVector(1.0, 1.0^2)
v3 = SVector(0.50, 0.50^2)

verts_parabola = [v1, v2, v3]
face_parabola = [SVector(1,2,3)]

m_parabola = CompScienceMeshes.CurvilinearMesh(verts_parabola, face_parabola, 1, 2)
ch_parabola = chart(m_parabola, 1)

mp_1 = neighborhood(ch_parabola, 0.0)
mp_2 = neighborhood(ch_parabola, 0.25)
mp_3 = neighborhood(ch_parabola, 0.50)
mp_4 = neighborhood(ch_parabola, 0.75)
mp_5 = neighborhood(ch_parabola, 1.0)

@test cartesian(mp_1) == SVector(1.00, 1.0)
@test cartesian(mp_2) == SVector(0.75, 0.5625)
@test cartesian(mp_3) == SVector(0.50, 0.25)
@test cartesian(mp_4) == SVector(0.25, 0.0625)
@test cartesian(mp_5) == SVector(0.00, 0.0)

@test tangents(mp_1, 1) == SVector(-1.0, -2.0)
@test tangents(mp_2, 1) ≈ SVector(-1.0, -1.5)
@test tangents(mp_3, 1) ≈ SVector(-1.0, -1.0)
@test tangents(mp_4, 1) ≈ SVector(-1.0, -0.5)
@test tangents(mp_5, 1) ≈ SVector(-1.0, -0.0)

@test volume(ch_parabola) ≈ 1.478942857544597 # (Reference from Mathematica 14)



## Triangle order 1
v1 = SVector(1.0, 0.0, 0.0)
v2 = SVector(0.0, 1.0, 0.0)
v3 = SVector(0.0, 0.0, 0.0)

verts1 = [v1, v2, v3]
faces1 = [SVector(1, 2, 3)]

m1 = CompScienceMeshes.Mesh(verts1, faces1)
m2 = CompScienceMeshes.CurvilinearMesh(verts1, faces1, 2, 1)

ch1 = chart(m1, 1)
ch2 = chart(m2, 1)

@test barytocart(ch1, SVector(1.0, 0.0)) == SVector(1.0, 0.0, 0.0)
@test barytocart(ch1, SVector(0.0, 1.0)) == SVector(0.0, 1.0, 0.0)
@test barytocart(ch1, SVector(0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
@test barytocart(ch2, SVector(1.0, 0.0)) == SVector(1.0, 0.0, 0.0)
@test barytocart(ch2, SVector(0.0, 1.0)) == SVector(0.0, 1.0, 0.0)
@test barytocart(ch2, SVector(0.0, 0.0)) == SVector(0.0, 0.0, 0.0)

CompScienceMeshes.unitarybase(ch2, SVector(0.0, 0.0))
CompScienceMeshes.unitarybase(ch2, SVector(1.0, 0.0))
CompScienceMeshes.unitarybase(ch2, SVector(0.0, 1.0))
CompScienceMeshes.unitarybase(ch2, SVector(0.5, 0.5))

## Verify consistency with linear mesh
v1 = SVector(0.0, 0.0)
v2 = SVector(1.0, 0.0)
v3 = SVector(0.25, 0.0)
v4 = SVector(0.50, 0.0)
v5 = SVector(0.75, 0.0)

verts1 = [v1, v2]
faces1 = [SVector(1,2)]

verts2 = [v1, v2, v3, v4, v5]
faces2 = [SVector(1,2,3,4,5)]

m1 = CompScienceMeshes.Mesh(verts1, faces1)
m2 = CompScienceMeshes.CurvilinearMesh(verts2, faces2, 1, 4)

ch1 = chart(m1, 1)
ch2 = chart(m2, 1)

@test barytocart(ch1, 1.0) == SVector(0.0, 0.0)
@test barytocart(ch2, 1.0) == SVector(0.0, 0.0)
@test barytocart(ch1, 0.0) == SVector(1.0, 0.0)
@test barytocart(ch2, 0.0) == SVector(1.0, 0.0)

mp1 = neighborhood(ch1, 0.0)
mp2 = neighborhood(ch2, 0.0)

@test cartesian(mp1) == SVector(1.0, 0.0)
@test cartesian(mp2) == SVector(1.0, 0.0) # We are using barycentric coordinates: 0.0 refers to the end vertex

@test tangents(mp1, 1) == SVector(-1.0, 0.0)
@test tangents(mp2, 1) == SVector(-1.0, 0.0) # Sign is flipped

@test normal(mp1) == SVector(0.0, -1.0)
@test normal(mp2) == SVector(0.0, -1.0) # Sign is flipped

@test nodes(ch1)[1] == SVector(0.0, 0.0)
@test nodes(ch1)[2] == SVector(1.0, 0.0)

@test nodes(ch2)[1] == SVector(0.0, 0.0)
@test nodes(ch2)[2] == SVector(1.0, 0.0)

@test volume(ch2) ≈ 1.0

##
v1 = SVector(1.0, 0.0, 0.0)
v2 = SVector(0.0, 1.0, 0.0)
v3 = SVector(0.0, 0.0, 0.0)

verts1 = [v1, v2, v3]
faces1 = [SVector(1, 2, 3)]


verts2 = [
    SVector(1.0, 0.0, 0.0),                         # 1
    SVector(0.0, 1.0, 0.0),                         # 2
    SVector(0.0, 0.0, 0.0),                         # 3
    SVector(0.0, 0.5000000000013304, 0.0),          # 4
    SVector(0.4999999999986718, 0.0, 0.0),          # 5
    SVector(0.5000000000013293, 0.4999999999986707, 0.0) # 6
]

face2 = [SVector(1, 2, 3, 6, 4, 5)]

m1 = CompScienceMeshes.Mesh(verts1, faces1)
m2 = CompScienceMeshes.CurvilinearMesh(verts2, face2, 2, 2)

ch1 = chart(m1, 1)
ch2 = chart(m2, 1)

## Order 7

v1 = SVector(1.0, 0.0, 0.0)
v2 = SVector(0.0, 1.0, 0.0)
v3 = SVector(0.0, 0.0, 0.0)

verts1 = [v1, v2, v3]
faces1 = [SVector(1, 2, 3)]

verts2 = [
    SVector(1.0, 0.0, 0.0),
    SVector(0.0, 1.0, 0.0),
    SVector(0.0, 0.0, 0.0),
    SVector(0.0, 0.8571428571429518, 0.0),
    SVector(0.0, 0.7142857142863779, 0.0),
    SVector(0.0, 0.5714285714298644, 0.0),
    SVector(0.0, 0.428571428572595, 0.0),
    SVector(0.0, 0.2857142857150634, 0.0),
    SVector(0.0, 0.1428571428575318, 0.0),
    SVector(0.1428571428568644, 0.0, 0.0),
    SVector(0.2857142857135911, 0.0, 0.0),
    SVector(0.4285714285703177, 0.0, 0.0),
    SVector(0.5714285714274409, 0.0, 0.0),
    SVector(0.7142857142849605, 0.0, 0.0),
    SVector(0.8571428571424802, 0.0, 0.0),
    SVector(0.8571428571430437, 0.1428571428569563, 0.0),
    SVector(0.7142857142863934, 0.2857142857136066, 0.0),
    SVector(0.5714285714297733, 0.4285714285702267, 0.0),
    SVector(0.428571428572577, 0.571428571427423, 0.0),
    SVector(0.2857142857150514, 0.7142857142849486, 0.0),
    SVector(0.1428571428575258, 0.8571428571424742, 0.0),
    SVector(0.7142857142855676, 0.1428571428569817, 0.0),
    SVector(0.1428571428574321, 0.7142857142854938, 0.0),
    SVector(0.1428571428569261, 0.1428571428574504, 0.0),
    SVector(0.5714285714289259, 0.2857142857136924, 0.0),
    SVector(0.4285714285721156, 0.4285714285706285, 0.0),
    SVector(0.2857142857148444, 0.5714285714280445, 0.0),
    SVector(0.1428571428573385, 0.5714285714289371, 0.0),
    SVector(0.1428571428571978, 0.4285714285721468, 0.0),
    SVector(0.1428571428570315, 0.2857142857148569, 0.0),
    SVector(0.285714285713734, 0.1428571428573691, 0.0),
    SVector(0.428571428570712, 0.1428571428572502, 0.0),
    SVector(0.571428571428096, 0.1428571428570705, 0.0),
    SVector(0.4285714285714073, 0.2857142857140783, 0.0),
    SVector(0.2857142857145434, 0.4285714285713746, 0.0),
    SVector(0.2857142857140778, 0.2857142857145756, 0.0)
]

face2 = [SVector(1, 2, 3, 16, 17, 18, 19, 20, 21, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36)]

m1 = CompScienceMeshes.Mesh(verts1, faces1)
m2 = CompScienceMeshes.CurvilinearMesh(verts2, face2, 2, 7)

ch1 = chart(m1, 1)
ch2 = chart(m2, 1)

@test barytocart(ch1, SVector(1.0, 0.0)) == SVector(1.0, 0.0, 0.0)
@test barytocart(ch1, SVector(0.0, 1.0)) == SVector(0.0, 1.0, 0.0)
@test barytocart(ch1, SVector(0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
@test barytocart(ch2, SVector(1.0, 0.0)) == SVector(1.0, 0.0, 0.0)
@test barytocart(ch2, SVector(0.0, 1.0)) == SVector(0.0, 1.0, 0.0)
@test barytocart(ch2, SVector(0.0, 0.0)) == SVector(0.0, 0.0, 0.0)

mp1 = neighborhood(ch1, SVector(0.0, 0.0))
mp2 = neighborhood(ch2, SVector(0.0, 0.0))

@test volume(ch1) == 0.5
@test volume(ch2) ≈ 0.5

@test cartesian(mp1) == SVector(0.0, 0.0, 0.0)
@test cartesian(mp2) == SVector(0.0, 0.0, 0.0) 

@test tangents(mp1, 1) == SVector(1.0, 0.0, 0.0)
@test tangents(mp1, 2) == SVector(0.0, 1.0, 0.0)

CompScienceMeshes.unitarybase(ch2, SVector(0.0, 0.0))
CompScienceMeshes.unitarybase(ch2, SVector(1.0, 0.0))
CompScienceMeshes.unitarybase(ch2, SVector(0.0, 1.0))

@test normal(ch2, SVector(0.0, 0.0)) == SVector(0.0, 0.0, 1.0)

@test tangents(mp2, 1) ≈ SVector(1.0, 0.0, 0.0) # Sign is flipped
@test tangents(mp2, 2) ≈ SVector(0.0, 1.0, 0.0) # Sign is flipped

@test normal(mp1) == SVector(0.0, 0.0, 1.0)
@test normal(mp2) == SVector(0.0, 0.0, 1.0)

@test jacobian(mp1) == 1.0
@test jacobian(mp2) ≈ 1.0
##

# Consider parabola t ↦(t, t²)
v1 = SVector(0.0, 0.0^2)
v2 = SVector(1.0, 1.0^2)
v3 = SVector(0.25, 0.25^2)
v4 = SVector(0.50, 0.50^2)
v5 = SVector(0.75, 0.75^2)

verts_parabola = [v1, v2, v3, v4, v5]
face_parabola = [SVector(1,2,3,4,5)]

m_parabola = CompScienceMeshes.CurvilinearMesh(verts_parabola, face_parabola, 1, 4)
ch_parabola = chart(m_parabola, 1)

mp_1 = neighborhood(ch_parabola, 0.0)
mp_2 = neighborhood(ch_parabola, 0.25)
mp_3 = neighborhood(ch_parabola, 0.50)
mp_4 = neighborhood(ch_parabola, 0.75)
mp_5 = neighborhood(ch_parabola, 1.0)

@test cartesian(mp_1) == SVector(1.00, 1.0)
@test cartesian(mp_2) == SVector(0.75, 0.5625)
@test cartesian(mp_3) == SVector(0.50, 0.25)
@test cartesian(mp_4) == SVector(0.25, 0.0625)
@test cartesian(mp_5) == SVector(0.00, 0.0)

@test tangents(mp_1, 1) == SVector(-1.0, -2.0) # Sign is flipped
@test tangents(mp_2, 1) ≈ SVector(-1.0, -1.5) # Sign is flipped
@test tangents(mp_3, 1) ≈ SVector(-1.0, -1.0) # Sign is flipped
@test tangents(mp_4, 1) ≈ SVector(-1.0, -0.5) # Sign is flipped
@test tangents(mp_5, 1) ≈ SVector(-1.0, -0.0) # Sign is flipped

@test volume(ch_parabola) ≈ 1.478942857544597 # (Reference from Mathematica 14)

#=
##
using CompScienceMeshes
using Test
using StaticArrays
fn = joinpath(dirname(@__FILE__),"triangle.msh")
m = CompScienceMeshes.load_gmsh_mesh(fn,  order=5)

m.faces

# Test that we are not missing a node
for p = 1:10
    triplets = CompScienceMeshes.gmsh_triangle_index_to_triplet[p]
    @test length(triplets) == div((p+1)*(p+2), 2)
    #map = gmsh_triangle_index_to_triplet(p)
    #for (i, t) in enumerate(triplets)
    #    @test map[t] == i - 1  # Gmsh uses 0-based indexing
    #end
end

function gmsh_triplet_to_index(p::Int)
    triplets = gmsh_index_to_triplet(p)
    map = Dict{Tuple{Int,Int,Int}, Int}()
    for (i, t) in enumerate(triplets)
        map[t] = i - 1  # Gmsh uses 0-based indexing
    end
    return map
end

=#