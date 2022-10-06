using CompScienceMeshes

# G = readmesh("examples/sphere.in")
# edges = skeleton(G,1)
# faces = skeleton(G,2)
# connectMat = connectivity(edges,faces)
# NV, NE, NF = numvertices(G), numcells(edges), numcells(faces)
# nv = NV + NE
# verts = Array{vertextype(G)}(nv)
# G1 = bisecting_refinement(G)
# G2 = Loop_subdivision(G)
# G3 = Loop_subdivision(G2)
# G4 = Loop_subdivision(G3)
# writemesh(G2,"examples/sphere_subd1.in")
# writemesh(G3,"examples/sphere_subd2.in")
# writemesh(G4,"examples/sphere_subd3.in")
fn = joinpath(dirname(@__FILE__),"sphere_subd3.in")
#G = readmesh("/Users/Benjamin/Documents/sphere.in")
G = readmesh(fn)
# G = Loop_subdivision(G)
subdG=GSubdMesh(G)
nelem = length(subdG.elements)
area = 0.0
ng = 4
gpt,wt = gauss_points(ng)
pts = []
for e = 1:1
    C = chart(subdG,e)
    for ig = 1 : ng*ng
        gpta = gpt[ig,1]
        gptb = gpt[ig,2]
        w = wt[ig]
        gpxi = (1+gpta)*(1-gptb)/4.0
        gpeta = (1+gptb)/2.0
        w *= (1-gpeta)/4.0
        g = neighborhood(C,[gpxi,gpeta])
        jac = jacobian(g)
        global area += jac * w
        push!(pts, cartesian(g))
    end
end
print(" area = ")
print(area)

