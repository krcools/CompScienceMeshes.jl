import Base.isless

#export sortneighbors

function isless{N,T}(p::SVector{N,T}, q::SVector{N,T})
    for i = 1 : N
        if isless(p[i], q[i])
            return true
        end
    end
    return false
end
