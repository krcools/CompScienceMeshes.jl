"""
N: universe dimension
D: manifold dimension
T: coordinate type
"""
immutable RectangleChart{N,D,T}
  origin::SVector{N,T}
  tangents::SVector{D,SVector{N,T}}
  normal::SVector{N,T}
end

function RectangleChart(o, tgs)
  n = tgs[1] × tgs[2]
  n = normalize(n)
  RectangleChart(o, tgs, n)
end

immutable RectangleNBD{N,D,T}
  chart::RectangleChart{N,D,T}
  parm::SVector{D,T}
  cart::SVector{N,T}
end

cartesian(nbd::RectangleNBD) = nbd.cart

function neighborhood(rect::RectangleChart{3,2}, uv)
  o = rect.origin
  t1, t2 = rect.tangents
  u1, u2 = uv
  cart = o + u1*t1 + u2*t2
  RectangleNBD(rect, SVector(uv), cart)
end

getindex(nbd::RectangleNBD, i) = nbd.cart[i]
normal(nbd::RectangleNBD) = nbd.chart.normal

function quadpoints(chart::RectangleChart{3,2}, rule)

  T = eltype(chart.origin)

  o = chart.origin
  t1 = chart.tangents[1]
  t2 = chart.tangents[2]
  j = norm(t1 × t2)
  Q, P = rule
  U1, W1 = legendre(Q, zero(T), one(T))
  U2, W2 = legendre(P, zero(T), one(T))

  #pts = [o + u1 * t1 + u2 * t2 for u1 in U1 for u2 in U2]
  pts = [neighborhood(chart, (u1,u2)) for u1 in U1 for u2 in U2]
  wts = [j*w1*w2 for w1 in W1 for w2 in W2]

  collect(zip(pts,wts))
end
