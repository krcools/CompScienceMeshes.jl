import Base.isless

#export sortneighbors

function isless{N,T}(p::Vec{N,T}, q::Vec{N,T})
    for i = 1 : N
        if isless(p[i], q[i])
            return true
        end
    end
    return false
end
