using Test
using CompScienceMeshes
using StaticArrays

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
m2 = CompScienceMeshes.CurvilinearMesh(verts2, faces2, 4)

ch1 = chart(m1, 1)
ch2 = chart(m2, 1)

mp1 = neighborhood(ch1, 0.0)
mp2 = neighborhood(ch2, 0.0)

@test tangents(mp1, 1) == SVector(-1.0, 0.0)
@test tangents(mp2, 1) == SVector(-1.0, 0.0) # Sign is flipped

@test normal(mp1) == SVector(0.0, -1.0)
@test normal(mp2) == SVector(0.0, -1.0) # Sign is flipped

@test cartesian(mp1) == SVector(1.0, 0.0)
@test cartesian(mp2) == SVector(1.0, 0.0) # We are using barycentric coordinates: 0.0 refers to the end vertex

@test nodes(ch1)[1] == SVector(0.0, 0.0)
@test nodes(ch1)[2] == SVector(1.0, 0.0)

@test nodes(ch2)[1] == SVector(0.0, 0.0)
@test nodes(ch2)[2] == SVector(1.0, 0.0)

@test volume(ch2) ≈ 1.0


##

# Consider parabola t ↦(t, t²)
v1 = SVector(0.0, 0.0^2)
v2 = SVector(1.0, 1.0^2)
v3 = SVector(0.25, 0.25^2)
v4 = SVector(0.50, 0.50^2)
v5 = SVector(0.75, 0.75^2)

verts_parabola = [v1, v2, v3, v4, v5]
face_parabola = [SVector(1,2,3,4,5)]

m_parabola = CompScienceMeshes.CurvilinearMesh(verts_parabola, face_parabola, 4)
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