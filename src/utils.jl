import Base.isless
import Base.getindex
import Base.setindex!

function isless{N,T}(p::Vec{N,T}, q::Vec{N,T})
    for i = 1 : N
        if isless(p[i], q[i])
            return true
        end
    end
    return false
end

function getindex{T,N}(A::AbstractArray{T}, I::Vec{N})
    B = Array(T,N)
    for i in 1:N
        B[i] = A[I[i]]
    end
    return B
end

# getindex(A::AbstractArray, I::Vec{1}) = Vec(A[I[1]])
# getindex(A::AbstractArray, I::Vec{2}) = Vec(A[I[1]], A[I[2]])
# getindex(A::AbstractArray, I::Vec{3}) = Vec(A[I[1]], A[I[2]], A[I[3]])

function setindex!{T,N}(A::AbstractArray{T}, v, I::Vec{N})
    for i in 1:N
        A[I[i]] = v
    end
end

function setindex!{N}(A::AbstractArray, V::AbstractArray, I::Vec{N})
    for i in 1:N
        A[I[i]] = V[i]
    end
end

function searchsortedfirst(v, x, lo, hi, by, lt)
    lo = lo-1
    hi = hi+1
    @inbounds while lo < hi-1
        m = (lo+hi)>>>1
        if lt(by(v[m]), x)
        #if lt(o, v[m], x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

function colsearchsortedfirst(A, col)

    n = size(A,2)
    v =  collect(1:n)
    searchsortedfirst(v, col, 1, n, i->A[:,i], lexless)

end
