using ClusterTrees

"""
    Sort objects to lie on a space filling curve
"""
function sort_sfc(points)

    # A non-finite coordinate makes the box size NaN, so the tree would subdivide forever
    all(p -> all(isfinite, p), points) ||
        throw(ArgumentError("sort_sfc requires finite point coordinates."))

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

@testitem "sort_sfc" begin
    using StaticArrays
    pts = [rand(SVector{3,Float64}) for _ in 1:100]
    @test sort(CompScienceMeshes.sort_sfc(pts)) == 1:100
    pts[17] = SVector(NaN, 0.0, 0.0)
    @test_throws ArgumentError CompScienceMeshes.sort_sfc(pts)
end
