export isinside
export isinclosure

function isinside(cell, point)

  u = carttobary(cell, point)
  T = eltype(u)

  tol = eps(T) * 1e3
  w = one(T)
  for i in 1:dimension(cell)
    0+tol < u[i] < 1-tol || return false
    w -= u[i]
  end

  0+tol < w < 1-tol || return false
  return true

end

function isinclosure(cell, point)

  u = carttobary(cell, point)
  T = eltype(u)

  tol = eps(T) * 1e3
  w = one(T)
  for i in 1:dimension(cell)
    0-tol < u[i] < 1+tol || return false
    w -= u[i]
  end

  0-tol < w < 1+tol || return false
  return true

end
