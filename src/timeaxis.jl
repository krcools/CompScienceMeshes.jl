"""
This type conforms to the mesh interface but is specialised for the case of a segment of the real axis subdived in equally sized intervals. Typical use is a discretisation of the time axis.
"""
mutable struct SegmentedAxis{T}
    timestep::T
    numsteps::Int
end

numcells(ax::SegmentedAxis) = ax.numsteps
#cellvertices(ax::SegmentedAxis, i) = SVector(point((i-1)*ax.timestep), point(i*ax.timestep))

#Base.size(a::SegmentedAxis) = (a.numsteps,)
#Base.linearindexing(a::SegmentedAxis) = Base.LinearFast()
#Base.getindex(a::SegmentedAxis, i) = ((i-1)*a.timestep, i*a.timestep)

struct edgevertexgeo{T,P}

    mindist::T #distance between the vertex and enter line support of the edge
    dmin :: T #from the actual edge
    dmax :: T 
    tangent::P #direction along the edge
    originpoint::P #ortogonal projection of the vertex on the line support of the edge
    extint0::NTuple{2,T} #coordinates of the extreme of the original edge
end


function edgevertexinteraction(rp,r1,r2) #convention: positive direction from r2 to r1,
    r12=r1-r2
    l12=norm(r12)
    x=r12/l12 
    rp1=rp-r1
    #lp1=norm(rp1)
    b=-dot(rp1,x)
    r01=-b*x
    r0=r01+r1
    #rp0=rp-r0
    a=dot(r2-r0,x)

    L=norm(rp1-r01) 

    T = eltype(rp)
    m = norm(rp-r2)
    M = m
    z=zero(eltype(rp))
    if a*b <= z
       m=L
       abs(a)<abs(b) && (M=norm(rp1))
        
    else
        for q in (r2,r1)
            d = norm(rp-q)
            d < m && (m=d)
            d > M && (M=d)
        end
    end
    
    return edgevertexgeo(
                L,
                m,
                M,
                x,
                r0,
                (a, b),
                )
end

function anglecontribution(ξ,n,geo)#this is used to check if the projection of a vertex onto planar support of a triangle is inside the triangle
    edgevertgeo=geo

    a=edgevertgeo.extint0[1]
    b=edgevertgeo.extint0[2]
    r0=edgevertgeo.originpoint

    t = edgevertgeo.tangent
    mseg = cross(t, n)
    p = dot(r0-ξ,mseg)

    if abs(p)<eps(eltype(ξ))
        α=zero(eltype(ξ))
    else
        α=atan(b/p)-atan(a/p)
    end
    return α
end

function mindist2D_edg_edg(edg1,edg2)
    a1,a2=edg1[1],edg1[2]
    b1,b2=edg2[1],edg2[2]
    vertices=[a1,a2,b1,b2]
    T=eltype(a1)
    geo1=edgevertexinteraction(a1,b1,b2) 
    geo2=edgevertexinteraction(a2,b1,b2)
    geo3=edgevertexinteraction(b1,a1,a2)
    geo4=edgevertexinteraction(b2,a1,a2)
    geo=[geo1,geo2,geo3,geo4]

    x=geo[3].tangent
    xb=geo[1].tangent

    hdir=cross(xb,x)

    angletot=zero(eltype(a1))
    dminv=Vector{eltype(a1)}(undef, 4)
    #dmaxv=Vector{eltype(a1)}(undef, 4)
    ξ=Vector{typeof(a1)}(undef, 4)
    for j in 1:4
        dminv[j]=geo[j].dmin
        #dmaxv[j]=geo[j].dmax 
    end

    if norm(hdir) < (100)*eps(eltype(a1))
        dmin=min(dminv[1],dminv[2],dminv[3],dminv[4])
        
        #dmax=max(dmaxv[1],dmaxv[2],dmaxv[3],dmaxv[4])
        return dmin
    end

    n=hdir/norm(hdir)
    sgnn=[+1,-1,-1,+1]
        h=dot(a2-b2,n)
        sgnh=[+1,-1,+1,-1]
        for j in 1:4 
            v=vertices[j]
            ξ[j]=v-n*h*sgnh[j]*sgnn[j] 
            angletot+=anglecontribution(ξ[j],sgnn[j]*n,geo[j])
        end
        if abs(angletot-2*T(π))<100*eps(eltype(a1))
            dmin=abs(h)
        else
            dmin=min(dminv[1],dminv[2],dminv[3],dminv[4])
        end

     #   dmax=max(dmaxv[1],dmaxv[2],dmaxv[3],dmaxv[4])
        return dmin

end

function mindist2D_nd_fc(node,face,geo)

    T=eltype(node)
    x=geo[3].tangent
    y=geo[2].tangent
    
    n=cross(x,y) 
    n/=norm(n)

    h=dot(node-face[1],n)
    ξ=node-n*h
    dminv=Vector{eltype(node)}(undef, 3)
    angletot=zero(eltype(node))
    for j in 1:3

        dminv[j]=geo[j].dmin
        #dmaxv[j]=geo[j].dmax 
            
        angletot+=anglecontribution(ξ,n,geo[j])
    end
    if abs(angletot-2*T(π))<100*eps(eltype(node))
        dmin=abs(h)
    else 
        dmin=min(dminv[1],dminv[2],dminv[3])
    end

    #dmax=max(dmaxv[1],dmaxv[2],dmaxv[3])
    return dmin
end    



function minmaxdist(τ, σ)
	
    T = eltype(τ[1])

    # Extract the vertices of the triangles τ and σ as tuples
    τ_vertices = (τ[1], τ[2], τ[3])
    σ_vertices = (σ[1], σ[2], σ[3])

    # Extract the edges of the triangles τ and σ as tuples
    τ_edges = (simplex(τ[1], τ[2]), simplex(τ[2], τ[3]), simplex(τ[3], τ[1]))
    σ_edges = (simplex(σ[1], σ[2]), simplex(σ[2], σ[3]), simplex(σ[3], σ[1]))

    # Compute the rings for each vertex of τ with respect to σ
    geo_τ_σ = ntuple(i -> ntuple(j -> edgevertexinteraction(τ_vertices[i], σ_edges[j][1], σ_edges[j][2]), 3), 3)
    geo_σ_τ = ntuple(i -> ntuple(j -> edgevertexinteraction(σ_vertices[i], τ_edges[j][1], τ_edges[j][2]), 3), 3)

    # Compute min distances node-to-face
    mindist_nodeτ_faceσ = ntuple(i -> mindist2D_nd_fc(τ_vertices[i], σ, (geo_τ_σ[i][1], geo_τ_σ[i][2], geo_τ_σ[i][3])), 3)
    mindist_nodeσ_faceτ = ntuple(i -> mindist2D_nd_fc(σ_vertices[i], τ, (geo_σ_τ[i][1], geo_σ_τ[i][2], geo_σ_τ[i][3])), 3)

    #compute min distances edge-to-edge
    mindist_edgτ_edgσ_vector = ntuple(k -> begin
    i = fld1(k, 3)  # row index (1-based)
    j = (k - 1) % 3 + 1  # column index (1-based)
    mindist2D_edg_edg(τ_edges[i], σ_edges[j])
        end, 9)


    # Compute the overall minimum
    min_dist = minimum((mindist_nodeτ_faceσ..., mindist_nodeσ_faceτ..., mindist_edgτ_edgσ_vector...))

    # Compute the maximum vertex-to-vertex distance
    max_dist = zero(T)
    for p in τ_vertices
        for q in σ_vertices
            dist = norm(p - q)
            if dist > max_dist
                max_dist = dist
            end
        end
    end

    return min_dist, max_dist
end

function rings(τ, σ, ΔR)
	m, M = minmaxdist(τ, σ)
	r0 = floor(Int, m/ΔR) + 1
	r1 = ceil(Int, M/ΔR)
	r0 : r1
end

ring(r, ΔR) = ((r-1)*ΔR, r*ΔR)
