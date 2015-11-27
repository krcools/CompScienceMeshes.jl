CSM = CompScienceMeshes

mesh = CSM.meshrectangle(1.0, 1.0, 0.5)
pred(simplex) = maximum(abs(mesh.vertices[1,simplex])) < eps(Float64)
sm = CSM.submesh(pred, mesh, 1)
@test sm.faces == [1 2; 2 3]
