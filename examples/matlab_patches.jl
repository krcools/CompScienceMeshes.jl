using MATLAB

function patch(Γ)
    # V = Float64[v[j] for v in Γ.vertices, j in 1:universedimension(Γ)]
    # F = Int[f[j] for f in Γ.faces, j in 1:dimension(Γ)+1]
    V = vertexarray(Γ)
    F = cellarray(Γ)

    @mput V F
    @matlab begin
        #M = Mesh(V,F)
        hold("on")
        #patch(M)
        patch("Vertices", V, "Faces", F, "FaceColor", "red")
    end

    mat"axis equal"
    mat"colorbar"
    mat"view(3)"

end

function patch(Γ, C)

    V = vertexarray(Γ)
    F = cellarray(Γ)

    @mput V F C
    @matlab begin
        #M = Mesh(V,F)
        hold("on")
        #patch(M, C, "LineStyle", "none")
        patch("Vertices", V, "Faces", F, "FaceVertexCData", C, "FaceColor", "flat")
    end

    mat"axis equal"
    mat"colorbar"
    mat"view(3)"
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
    @matlab quiver3(C(:,1),C(:,2),C(:,3), N(:,1),N(:,2),N(:,3))
end

function jmatlab_plot(x,y)
    @mput x
    @mput y
    @matlab plot(x,y)
end


# function quiver(Γ, fcr)
#
#   D = length(eltype(fcr))
#   Q = Float64[real(f[i]) for f in fcr, i in 1:D]
#   V = vertexarray(Γ)
#   F = cellarray(Γ)
#
#   @mput V F Q
#   mat"""
#   M = Mesh(V,F);
#   [G, N] = faceNormals(M);
#   quiver3x(G,Q,3);
#   """
# end
