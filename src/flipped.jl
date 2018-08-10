struct FlippedMesh{M}
        base::M
end

Base.:-(m::FlippedMesh) = m.base

struct FlippedMeshCellIterator{CA}
    cellarray::CA
end

cells(m::FlippedMesh) = FlippedMeshCellIterator(m.base.faces)

function Base.iterate(iter::FlippedMeshCellIterator)
    it = iterate(iter.cellarray)
    (it == nothing) ? it : (flip(it[1]), it[2],)
end

function Base.iterate(iter::FlippedMeshCellIterator, state)
    it = iterate(iter.cellarray, state)
    (it == nothing) ? it : (flip(it[1]), it[2],)
end

# Base.start(iter::FlippedMeshCellIterator) = start(iter.cellarray)
# Base.done(iter::FlippedMeshCellIterator, state) = done(iter.cellarray, state)
# function Base.next(iter::FlippedMeshCellIterator, state)
#     item, newstate = next(iter.cellarray, state)
#     flip(item), newstate
# end
