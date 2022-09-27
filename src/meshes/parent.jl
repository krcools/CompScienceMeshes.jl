struct BarcentricRefinedParent{U,D1,T} <: AbstractMesh{U,D1,T}
    parent::AbstractMesh{U,D1,T}
    refined::AbstractMesh{U,D1,T}
end

function children(mesh::BarcentricRefinedParent, cell_idx)
    dim = dimension(mesh)
    num_children = factorial(dim+1)
    fine_cells = cells(mesh.refined)
    [fine_cells[num_children*(cell_idx-1)+i] for i in 1:num_children]
end