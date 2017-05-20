abstract AbstractMapper

mapper(mesh) = Dict((c,i) for (i,c) in enumerate(cells(mesh)))
