import Base.isless
import Base.getindex
import Base.setindex!

export sortneighbors

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

"""
Move the s-th element right after the d-th
"""
function move_after!(p,n,s,d)

    n[d] == s && return

    t1 = n[d]
    t2 = n[s]
    t3 = p[s]

    n[d] = s
    n[s] = t1
    p[s] == 0 || (n[p[s]] = t2)

    p[s] = d
    t1 == 0 || (p[t1] = s)
    t2 == 0 || (p[t2] = t3)
end

function move_before!(p, n, s, d)
    p[d] == s && return

    t1 = p[d]
    t2 = p[s]
    t3 = n[s]

    p[d] = s
    p[s] = t1
    n[s] == 0 || (p[n[s]] = t2)

    n[s] = d
    t1 == 0 || (n[t1] = s)
    t2 == 0 || (n[t2] = t3)
end

# function splice(p, n, s0, s1, d)
#     p[d] == s1 && return
#
#     t1 = p[d]
#     t2 = p[s0]
#     t3 = n[s1]
#
#     p[d] = s
#     p[s0] = t1
#     n[s1] == 0 || (p[n[s1]] = t2)
#
#     n[s1] = d
#     t1 == 0 || (n[t1] = s)
#     t2 == 0 || (n[t2] = t3)
# end

function sortneighbors(a, pred)

    n = collect(2:length(a)+1); n[end] = 0
    p = collect(0:length(a)-1); p[1] = 0

    last = 1
    while last != 0
        cand = n[last]
        while cand != 0
            pred(a[last], a[cand]) && break
            cand = n[cand]
        end
        cand == 0 && break
        move_after!(p, n, cand, last)
        last = cand
    end

    first = 1
    while last != 0
        cand = n[last]
        while cand != 0
            pred(a[cand], a[first]) && break
            cand = n[cand]
        end
        cand == 0 && break
        move_before!(p, n, cand, first)
        first = cand
    end

    b = similar(a)
    i, j = last, length(n)
    while true
        b[j] = a[i]
        i = p[i]
        i == 0 && break
        j -= 1
    end

    return b

end

isclosed(a, pred) = pred(a[end], a[1])
