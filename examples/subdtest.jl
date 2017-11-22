using CompScienceMeshes
G = readmesh("/Users/Benjamin/Documents/sphere.in")
subdG=GSubdMesh(G)
nelem = length(subdG.elements)
area = 0.0
ng = 10
gpt,wt = gauss_points(ng)
for e = 1:nelem
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
        area += jac * w
    end
end
print(" area = ")
print(area)
