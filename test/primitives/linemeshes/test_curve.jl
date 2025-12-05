
using CompScienceMeshes
using StaticArrays
using Test

curve(t) = SVector((9/20-(1/9)cos(5t))*cos(t), (9/20-(1/9)cos(5t))*sin(t))

# Detect closed curve
msh = meshcurve(curve, 0.1; tend=Float64(2π))
@test msh.faces[end][2] == 1

# Detect open curve
msh = meshcurve(curve, 0.1; tend=Float64(π))
@test msh.faces[end][2] != 1