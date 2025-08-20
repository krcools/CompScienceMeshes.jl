test_minmaxdist
using Test
using CompScienceMeshes
using LinearAlgebra
m1= meshcuboid(1.0,1.0,1.0,0.4)
#m1=meshsphere(1.0,h=0.3)


ch_faces=[chart(m1,i) for i in 1:length(m1)]

function inside_shell(x,y,Rmin,Rmax)
    R=norm(x-y)
    if Rmin < R < Rmax
        return one(eltype(x))
    else 
        return zero(eltype(x))
    end
end 

function quadrule_inshell(face_test,face_trial,Rmin,Rmax)
    T = face_test
    S= face_trial
    
    #testing quadrature rules
       
    
    
        #testing quadrature rules
    
        qpsT = quadpoints(T,8)
        qpsS= quadpoints(S,8)
        
    
        r = zero(eltype(T[1]))
    for qpS in qpsS

        pS, wS= qpS
        for qpT in qpsT
    
            pT, wT = qpT
            x = cartesian(pT)
            y= cartesian(pS)
    
    
                r +=wT * wS *inside_shell(x,y,Rmin,Rmax)
    
    
        end
    end
        return r
    
end
Ty=eltype(ch_faces[1][1])


function old_minmax_dist(τ, σ)
	
    T = eltype(τ[1])

    τ_vertices = (τ[1], τ[2], τ[3])
    σ_vertices = (σ[1], σ[2], σ[3])

    max_dist = zero(T)
    min_dist = norm(τ[1]-σ[1])
    for p in τ_vertices
        for q in σ_vertices
            dist = norm(p - q)
            if dist > max_dist
                max_dist = dist
            end
            if dist < min_dist
                min_dist = dist
            end
        end
    end
    return min_dist,max_dist
end

for i in 1:length(ch_faces)
    for j in 1:length(ch_faces)
        dmin,dmax=minmaxdist(ch_faces[i],ch_faces[j])
        old_dmin,old_dmax=old_minmax_dist(ch_faces[i],ch_faces[j])
        @test old_dmax ≈ dmax
        if abs(dmin-old_dmin)>100*eps(Ty)
            @test quadrule_inshell(ch_faces[i],ch_faces[j],dmin/2-15*eps(Ty),dmin-10*eps(Ty)) ≈ zero(Ty)
            res=quadrule_inshell(ch_faces[i],ch_faces[j],dmin,old_dmin)
            #@test  res > 1000*eps(Ty)
            #@show dmin, old_dmin
            #@show res
            #@test quadrule_inshell(ch_faces[i],ch_faces[j],dmax+10*eps(Ty),dmax+dmax/2+10*eps(Ty)) == zero(Ty)
        end
    end
end