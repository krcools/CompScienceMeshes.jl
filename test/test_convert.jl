using Test
using CompScienceMeshes
using LinearAlgebra
using StaticArrays

Γ32 = meshrectangle(Float32(1.0), Float32(1.0), Float32(0.2))
Γ64 = meshrectangle(Float64(1.0), Float64(1.0), Float64(0.2))

Γ_converted = convert(Float32, Γ64)

@test Γ32.vertices == Γ_converted.vertices
@test Γ32.faces == Γ_converted.faces
