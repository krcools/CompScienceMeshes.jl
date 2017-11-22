using MATLAB

function patch(Γ)

    V = vertexarray(Γ)
    F = cellarray(Γ)

    @mput V F
    mat"""
        hold('on')
        patch('Vertices', V, 'Faces', F, 'FaceColor', 'red')
    """

    mat"axis equal"
    mat"colorbar"
    mat"view(3)"

end

function patch(Γ, C)

    V = vertexarray(Γ)
    F = cellarray(Γ)

    @mput V F C
    mat"""
        hold('on')
        patch('Vertices', V, 'Faces', F, 'FaceVertexCData', C, 'FaceColor', 'flat')
    """

    mat"axis equal"
    mat"colorbar"
    mat"view(3)"
end


function jmatlab_quiver(mesh, facecurs)
    els = elements(mesh)
    cts = center.(els)
    C = [ct[i] for ct in cts, i in 1:3]
    N = [fc[i] for fc in facecurs, i in 1:3]
    @mput C N
    mat"quiver3(C(:,1),C(:,2),C(:,3), N(:,1),N(:,2),N(:,3))"
end

function jmatlab_quiver(m)
    els = [chart(m, cells(m,i)) for i in 1:numcells(m)]
    C = zeros(numcells(m),3)
    N = zeros(numcells(m),3)
    for i in 1:size(C,1)
        s = els[i]
        p = neighborhood(s, [1,1]/3)
        n = normal(p)
        c = cartesian(p)
        C[i,:] = [c[1],c[2],c[3]]
        N[i,:] = [n[1],n[2],n[3]]
    end
    @mput C N
    mat"quiver3(C(:,1),C(:,2),C(:,3), N(:,1),N(:,2),N(:,3))"
end

function jmatlab_plot(x,y)
    @mput x
    @mput y
    mat"plot(x,y)"
end


# function quiver(Γ, fcr)
#
#   D = length(eltype(fcr))
#   V = vertexarray(Γ)
#   F = cellarray(Γ)
#
#   C = [cartesian(center(chart(Γ,c))) for c in cells(Γ)]
#   C = [C[i][j] for i in eachindex(C), j in 1:3]
#
#   Q = [fcr[i][j] for i in eachindex(fcr), j in 1:3]
#
#   @show size(C)
#   @show size(Q)
#
#   @mput V F C Q
#   mat"""
#   patch('Vertices', V, 'Faces', F, 'FaceColor', 'red')
#   quiver3(C(:,1),C(:,2),C(:,3),Q(:,1),Q(:,2),Q(:,3))
#   """
# end
