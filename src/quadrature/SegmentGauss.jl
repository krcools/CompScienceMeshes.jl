import FastGaussQuadrature:gausslegendre

export legendre

function legendre(n::Integer, a::T, b::T) where T<:AbstractFloat
  xw = gausslegendre(n)
  x = xw[1]
  w = xw[2]
  s = (b-a) * (x+1) / 2 + a
  v = (b-a) / 2 * w
  return Array(s'), v
end
