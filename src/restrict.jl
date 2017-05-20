"""
    restriction(submesh, supermesh)

Computes the restriction matrix relative to a submesh `submesh` of `supermesh`.

The restriction matrix has size `(m,n)`, where

    m == numcells(submesh)
    n == numcells(supermesh)

It has entries `1` at location `[i,j]` iff cell `i` of submesh equals cell `j` of supermesh.
The remaining entries are zero.

This matrix is referred to as the restriction matrix because if it acts on an array of
samples taken at the cells of `supermesh` is selects out the samples in the cells that
are retained in `submesh`, taking into account any renumbering. Likewise, its transpose
is sometimes referred to as the extension-by-zero operator because it maps arrays of samples
taken in the cells of `submesh` into an array of samples taken in the cells of `supermesh`
by inserting zeros at cells that were not retained in `submesh`.
"""
function restriction(submesh::Mesh, supermesh::Mesh)
    index_in_supermesh = mapper(supermesh)
    R = spzeros(Int, numcells(submesh), numcells(supermesh))
    for (i,c) in enumerate(cells(submesh))
        j = index_in_supermesh[c]
        R[i,j] = 1
    end
    R
end
