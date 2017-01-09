using MATLAB

function patch(Γ)
    # V = Float64[v[j] for v in Γ.vertices, j in 1:universedimension(Γ)]
    # F = Int[f[j] for f in Γ.faces, j in 1:dimension(Γ)+1]
    V = vertexarray(Γ)
    F = cellarray(Γ)

    @mput V F
    @matlab begin
        M = Mesh(V,F)
        hold("on")
        patch(M)
    end

    mat"axis equal"
    mat"colorbar"
    mat"view(3)"

end

function patch(Γ, fcr)
    C = fcr #Float64[sqrt(real(dot(f,f))) for f in fcr]
    # V = Float64[v[j] for v in Γ.vertices, j in 1:3]
    # F = Int[f[j] for f in Γ.faces, j in 1:3]
    V = vertexarray(Γ)
    F = cellarray(Γ)

    @mput V F C
    @show C
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

function quiver(Γ, fcr)

  D = length(eltype(fcr))
  Q = Float64[real(f[i]) for f in fcr, i in 1:D]
  V = vertexarray(Γ)
  F = cellarray(Γ)

  @mput V F Q
  mat"""
  M = Mesh(V,F);
  [G, N] = faceNormals(M);
  quiver3x(G,Q,3);
  """
end
