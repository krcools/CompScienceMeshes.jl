using CompScienceMeshes

"""
    Creates a predicate that checks that 2 edges belonging to `small` and `big`
    lie on the same plane and the normals of its surrounding cells points in the
    same direction
"""
function normal_cpredicate(small, big)
    @assert dimension(small) == 2
    @assert dimension(big) == 2

    ss = skeleton(small,1)
    sb = skeleton(big,1)

    D1 = connectivity(ss,small)
    D2 = connectivity(sb ,big)

    function pred(e1,e2)
        #TODO:Modify to take simplex rather than edge indices
        edgefaces1 = D1[:,e1]
        edgefaces2 = D2[:,e2]
        faceid1 = edgefaces1.nzind
        faceid2 = edgefaces2.nzind

        N1 = length(faceid1)
        N2 = length(faceid2)
        narray1 = Vector{CompScienceMeshes.StaticArrays.SArray{
                Tuple{3},Float64,1,3}}()
        narray2 = copy(narray1)

        for i = 1:N1
            push!(narray1,chart(small,small.faces[faceid1[i]]).normals[1])
        end
        @assert all(y->y==narray1[1],narray1) "egde on small mesh found to lie on junction"

        for i = 1:N2
            push!(narray2,chart(big,big.faces[faceid2[i]]).normals[1])
        end

        @assert all(y->y==narray2[1],narray2) "egde on big mesh found to lie on junction"

        return narray1[1] == narray2[1]
    end
end

"""
    Enforces extra constraint on edges in `small` and `big` that overlap with
    edges in `pInt` during `embedding` to ensure that both edges in `small` and
    `big` lie on planes with normals pointing in the same direction.
"""
function reembedding(small, big, pInt)
    @assert dimension(small) == 2
    @assert dimension(big) == 2
    @assert dimension(pInt) == 2

    ss = skeleton(small,1)
    sb = skeleton(big,1)
    sp = CompScienceMeshes.interior(pInt)

    epred = normal_cpredicate(small,big)
    overlaps = overlap_gpredicate(sp)
    mpred(e1,e2) = overlaps(chart(sb,cells(sb)[e2])) ? epred(e1,e2) : true
    #applies constraint during embedding only on edges that overlap with pInt
    CompScienceMeshes.embedding(ss,sb,predicate=mpred)
end


isclosed(m) = (numcells(boundary(m)) == 0)

plate = meshrectangle(1.0,1.0,1.0)
G01 = CompScienceMeshes.SComplex2D(plate)
PT = weld(G01, -G01, seam=boundary(G01));

A1 = reembedding( G01, PT, G01)
A2 = reembedding(-G01, PT, G01)

@show A1
@show A2
