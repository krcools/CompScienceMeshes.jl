using ClusterTrees

"""
    Sort objects to lie on a space filling curve
"""
function sort_sfc(points)

    ct, sz = CompScienceMeshes.boundingbox(points)
    tree = ClusterTrees.LevelledTrees.LevelledTree(ct, sz, Int[])

    smb = sz / 2^log(length(points)+1)
    for (i,pt) = enumerate(points)
        dest = (smallest_box_size=smb, target_point=pt)
        state = ClusterTrees.LevelledTrees.rootstate(tree, dest)
        ClusterTrees.update!(tree, state, i, dest) do tree, node, i
            push!(data(tree, node).values, i)
        end
    end

    sorted = Vector{Int}()
    for node in ClusterTrees.DepthFirstIterator(tree, root(tree))
        append!(sorted, data(tree,node).values)
    end

    return sorted
end
