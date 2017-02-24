"""
This type conforms to the mesh interface but is specialised for the case of a segment of the real axis subdived in equally sized intervals. Typical use is a discretisation of the time axis.
"""
type SegmentedAxis{T}
    timestep::T
    numsteps::Int
end

numcells(ax::SegmentedAxis) = ax.numsteps
cellvertices(ax::SegmentedAxis, i) = SVector(point((i-1)*ax.timestep), point(i*ax.timestep))

Base.size(a::SegmentedAxis) = (a.numsteps,)
Base.linearindexing(a::SegmentedAxis) = Base.LinearFast()
Base.getindex(a::SegmentedAxis, i) = ((i-1)*a.timestep, i*a.timestep)

function minmaxdist(τ, σ)
	T = eltype(τ[1])
	m = norm(τ[1]-σ[1])
	M = zero(T)
	for i in 1:length(τ)
		p = τ[i]
		for j in 1:length(σ)
			q = σ[j]
			d = norm(p-q)
			d < m && (m=d)
			d > M && (M=d)
		end
	end
	return m, M
end

function rings(τ, σ, ΔR)
	m, M = minmaxdist(τ, σ)
	r0 = floor(Int, m/ΔR) + 1
	r1 = ceil(Int, M/ΔR+1)
	r0 : r1
end

ring(r, ΔR) = ((r-1)*ΔR, r*ΔR)
