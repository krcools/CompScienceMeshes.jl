using CompScienceMeshes
using SauterSchwabQuadrature
#create a subdivision surfaces for a sphere exaple
G = readmesh("/Users/Benjamin/Documents/sphere.in")
subdG=GSubdMesh(G)

chart = chart(subdG,1)
u1= [0.5,0.25]
u2= [0.7,0.1]
d1 = neighborhood(chart,u1)
d2 = neighborhood(chart,u2)
global p1 = cartesian(d1)
global p2 = cartesian(d2)
function kernel(x,y)
			return(((x-p1)'*(y-p2))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end
result = sauterschwabintegral(chart, kernel)
