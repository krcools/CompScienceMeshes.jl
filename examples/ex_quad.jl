using CompScienceMeshes

o,x,y,z = CompScienceMeshes.euclidianbasis(3)

tet = simplex(o,x,2y,3z)
volume(tet)

qw = quadpoints(tet,4)
length(qw)

sum(w for (p,w) in qw) 
