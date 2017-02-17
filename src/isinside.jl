"""
    isinside(chart, point) -> Bool

Returns true is the given point is in the image of the given chart, false otherwise.
"""
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


"""
	isinclosure(simplex, point) -> Bool

Determine whether point is in the closure of simplex. False positives are possible for points just outside of the simplex.
"""
function isinclosure(cell, point)

  u = carttobary(cell, point)
  T = eltype(u)
  dim = dimension(cell)
  udim = length(point)

  tol = eps(T) * 1e3
  w = one(T)
  for i in 1:dim
    0-tol < u[i] < 1+tol || return false
    w -= u[i]
  end

  0-tol < w < 1+tol || return false

  # finally, check that the point indeed is in the plane of the cell
  A = zeros(T, udim, dim+1)
  origin = cell[dim+1]
  for i in 1:dim
      for j in 1:udim
          A[j,i] = cell[i][j] - origin[j]
      end
  end
  for j in 1:udim; A[j,end] = point[j] - origin[j]; end
  rank(A) == dim+1 && return false

  return true

end
