import Base.start
import Base.next
import Base.done

export normal
export cartesian, jacobian, unormal, meshpoint
export utangent, cellenumeration

export numcells
function numcells{M<:Mesh}(geo::Vector{M})
    nc = 0
    for mesh in geo
        nc += numcells(mesh)
    end
    return nc
end



function normal(vectors::AbstractArray)

    dim = size(vectors, 1)

    # we need an hyperplane as input
    @assert size(vectors,2) == dim-1

    n = Array(eltype(vectors), dim)
    for i in 1 : dim
        b = sub(vectors, [1:i-1; i+1:dim], :)
        n[i] = (-1)^(i-1) * det(b)
    end

    return (-1)^(dim+1) * n
end



immutable CellEnumerator{GeoType}
    geometry::GeoType
end

immutable CellEnumeratorState
    meshid::Int
    cellid::Int
    counter::Int
end

cellenumeration(geometry) = CellEnumerator(geometry)
start(cellenum::CellEnumerator) = CellEnumeratorState(1,1,1)

function next(cellenum::CellEnumerator, state)

    geometry = cellenum.geometry
    mesh  = geometry[state.meshid]

    index = state.counter
    indcs = mesh.faces[state.cellid]
    verts = mesh.vertices[indcs]
    #value = patch(verts, Val{CompScienceMeshes.dimension(mesh)})
    value = simplex(verts, Val{dimension(mesh)})

    newcellid = state.cellid + 1
    newmeshid = state.meshid
    if newcellid > numcells(mesh)
        newcellid = 1
        newmeshid += 1
    end

    newstate = CellEnumeratorState(newmeshid, newcellid, state.counter + 1)

    (index, value), newstate
end

function done(cellenum::CellEnumerator, state)
    geometry = cellenum.geometry
    if state.meshid > length(geometry) return true end
    if state.cellid > numcells(geometry[state.meshid]) return true end
    false
end

# not exported
