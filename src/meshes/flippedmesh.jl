struct FlippedMesh{U,D1,T} <: AbstractMesh{U,D1,T}
    mesh::AbstractMesh{U,D1,T}
end
function ==(m1::FlippedMesh{U,D1,T}, m2::FlippedMesh{U,D1,T}) where {U,D1,T}
    return m1.mesh == m2.mesh
end

indices(m::FlippedMesh, i) = fliporientation(indices(m.mesh,i))

cells(m::FlippedMesh) = (fliporientation(c) for c in cells(m.mesh))
vertices(m::FlippedMesh) = m.mesh.vertices

Base.:-(m::AbstractMesh) = FlippedMesh(m)
Base.:-(m::FlippedMesh) = m.mesh

Base.iterate(m::FlippedMesh, state=0) = iterate(m.mesh, state)
