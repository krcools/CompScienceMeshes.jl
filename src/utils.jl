import Base.isless

export sortneighbors

function isless{N,T}(p::Vec{N,T}, q::Vec{N,T})
    for i = 1 : N
        if isless(p[i], q[i])
            return true
        end
    end
    return false
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
