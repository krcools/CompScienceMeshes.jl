import PyPlot
import CompScienceMeshes
CM = CompScienceMeshes

function plot(mesh::CM.Mesh{3,3,Float64})
    X = [mesh.vertices[i][1] for i in 1 : CM.numvertices(mesh)]
    Y = [mesh.vertices[i][2] for i in 1 : CM.numvertices(mesh)]
    Z = [mesh.vertices[i][3] for i in 1 : CM.numvertices(mesh)]
    T = [mesh.faces[i][j] for i in 1 : CM.numcells(mesh), j in 1:3] - 1
    PyPlot.plot_trisurf(X,Y,Z,triangles=T)
end

function plot(mesh::CM.Mesh{3,2,Float64})
    X = [mesh.vertices[i][1] for i in 1 : CM.numvertices(mesh)]
    Y = [mesh.vertices[i][2] for i in 1 : CM.numvertices(mesh)]
    Z = [mesh.vertices[i][3] for i in 1 : CM.numvertices(mesh)]
    N = [1,2,2]
    T = [mesh.faces[i][N[j]] for i in 1 : CM.numcells(mesh), j in 1:3] - 1
    PyPlot.plot_trisurf(X,Y,Z,triangles=T)
end
