abstract type AbstractMapper end

mapper(mesh) = Dict((c,i) for (i,c) in enumerate(cells(mesh)))
