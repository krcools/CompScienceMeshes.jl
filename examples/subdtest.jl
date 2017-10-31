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
G = readmesh("/Users/Benjamin/Documents/sphere.in")
# G = Loop_subdivision(G)
subdG=GSubdMesh(G)
# writemesh(G2,"/Users/Benjamin/Documents/sphere2.in")
subd_elements = subdG.elements
vertices = G.vertices
nelem = length(subd_elements)
area = 0.0
ng = 10
# gpt,wt = gauss_quad_linear(10,-1,1)
gpt,wt = gauss_points(ng)
for (e, subdElem) in enumerate(subd_elements)
    Nodes = subdElem.RingNodes
    val = length(Nodes) - 6
    vertices_patch = zeros(3,length(Nodes))
    for i = 1:length(Nodes)
        vertices_patch[:,i] = vertices[Nodes[i]]'
    end
    for ig = 1 : ng*ng
        gpta = gpt[ig,1]
        gptb = gpt[ig,2]
        w = wt[ig]
        gpxi = (1+gpta)*(1-gptb)/4.0
        gpeta = (1+gptb)/2.0
        w *= (1-gpeta)/4.0
            # shape = shape_function(gpxi,gpeta,val)
        shape_der = shape_function_der(gpxi,gpeta,val)
        t1 = vertices_patch * shape_der[:,1]
        t2 = vertices_patch * shape_der[:,2]
        norm = [t1[2] * t2[3] - t1[3] * t2[2];
		        t1[3] * t2[1] - t1[1] * t2[3];
				t1[1] * t2[2] - t1[2] * t2[1]]
        jac = sqrt(norm[1]*norm[1] + norm[2]*norm[2] + norm[3]*norm[3])
        area += jac * w
    end
end
print(" area = ")
print(area)
