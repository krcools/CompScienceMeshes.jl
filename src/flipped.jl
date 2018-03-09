struct FlippedMesh{M}
        base::M
end

Base.:-(m::FlippedMesh) = m.base

#cells(m::FlippedMesh,i) = flip(cells(m.base,i))

struct FlippedMeshCellIterator{CA}
    cellarray::CA
end

cells(m::FlippedMesh) = FlippedMeshCellIterator(m.base.faces)

Base.start(iter::FlippedMeshCellIterator) = start(iter.cellarray)
Base.done(iter::FlippedMeshCellIterator, state) = done(iter.cellarray, state)
function Base.next(iter::FlippedMeshCellIterator, state)
    item, newstate = next(iter.cellarray, state)
    flip(item), newstate
end
