immutable FlippedMesh{M}
        base::M
end

Base.:-(m::FlippedMesh) = m.base

cells(m::FlippedMesh,i) = flip(cells(m.base,i))
