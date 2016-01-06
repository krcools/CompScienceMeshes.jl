using FixedSizeArrays
import Base.dot
import Base.*

@generated function *{T, M, N}(a::Mat{M, N, T}, b::Point{N,T})
    expr = [:(bilindot(row(a, $i), b.(1))) for i=1:M]
    :( Vec($(expr...)) )
end

@generated function dot{N,T,U}(a::Point{N,U}, b::Point{N,T})
  ex = :(a[1]*b[1])
  for i in 2:N
    ex = :($ex + a[$i]*b[$i])
  end
  return ex
end

@generated function *{T,N}(a::AbstractArray{T,2}, v::Vec{N,T})
    expr = Expr(:call, :+, [:(a[i,$j]*v[$j]) for j in 1:N]...)
    return quote
        @assert size(a,2) == length(v)
        r = Array{T}(size(a,1))
        for i in 1:size(a,1)
            r[i] = $expr
        end
        return r
    end
end

@generated function *{T,M,N}(a::Mat{M,N,T}, v::AbstractArray{T,1})
    expr = Expr(:call, :Vec, [Expr(:call, :+, [:(a[$m,$n]*v[$n]) for n in 1:N]...) for m in 1:M]...)
    return quote
        @assert N == length(v)
        $expr
    end
end
