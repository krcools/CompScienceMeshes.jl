struct FlippedMesh{U,D1,T} <: AbstractMesh{U,D1,T}
    mesh::AbstractMesh{U,D1,T}
end

cells(m::FlippedMesh) = (fliporientation(c) for c in cells(m.mesh))
vertices(m::FlippedMesh) = m.mesh.vertices

Base.:-(m::AbstractMesh) = FlippedMesh(m)
Base.:-(m::FlippedMesh) = m.mesh
