function meshsegment{T<:Real}(L::T, delta::T, udim=2)

  PT = Vec{udim, T}
  CT = Vec{2,Int}

  num_segments = ceil(Int, L/delta)
  actual_delta = L/num_segments
  x = collect(0:num_segments) * actual_delta

  vertices = zeros(PT, num_segments+1)
  for i in 1 : length(vertices)
    a = zeros(T, udim)
    a[1] = x[i]
    vertices[i] = PT(a)
  end

  faces = Array(CT, num_segments)
  for i in 1 : length(faces)
    faces[i] = Vec(i, i+1)
  end

  Mesh(vertices, faces)
end

function meshcircle{T<:Real}(radius::T, delta::T, udim=2)

  PT = Vec{udim,T}
  CT = Vec{2,Int}

    circumf = 2 * pi *radius
    num_segments = ceil(Int, circumf / delta)
    delta = circumf / num_segments
    dalpha = delta / radius
    alpha = collect(0 : num_segments-1) * dalpha

    vertices = Array(PT, num_segments)
    for i in 1 : num_segments
        a = zeros(T, udim)
        a[1] = radius * cos(alpha[i])
        a[2] = radius * sin(alpha[i])
        vertices[i] = PT(a)
    end

    faces = Array(CT, num_segments)
    for i in 1 : length(faces)-1
        faces[i] = Vec{2,Int}(i, i+1)
    end
    faces[end] = Vec{2,Int}(num_segments, 1)

    return Mesh(vertices, faces)
end

function meshrectangle{T}(width::T, height::T, delta::T, udim=3)

  PT = Vec{udim,T}
  CT = Vec{3,Int}

    @assert 2 <= udim

    #UDim = 3
    #MDim = 2

    nx = round(Int, ceil(width/delta));  nx = max(nx,1); dx = width/nx
    ny = round(Int, ceil(height/delta)); ny = max(ny,1); dy = height/ny

    xs = (0:nx) * dx
    ys = (0:ny) * dy

    vertices = zeros(PT, (nx+1)*(ny+1))
    k = 1
    for x in xs
        for y in ys
            p = zeros(T, udim)
            p[1] = x
            p[2] = y
            vertices[k] = PT(p)
            k += 1
        end
    end

    faces = zeros(CT, 2*nx*ny)
    k = 1
    for i in 1 : nx
        for j in 1 : ny
            v11 = (i-1)*(ny+1) + j
            v12 = i*(ny+1) + j
            v21 = (i-1)*(ny+1) + (j+1)
            v22 = i*(ny+1) + (j+1)
            faces[k]   = CT(v11,v21,v12)
            faces[k+1] = CT(v21,v22,v12)
            k += 2
        end
    end

    Mesh(vertices, faces)
end


function meshsphere(r, h)
    s =
"""
lc = $h;

Point(1)={0,0,0,lc};
Point(2)={$r,0,0,lc};
Point(3)={0,$r,0,lc};
Point(4)={-$r,0,0,lc};
Point(5)={0,-$r,0,lc};
Point(6)={0,0,$r,lc};
Point(7)={0,0,-$r,lc};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};
Circle(5)={6,1,3};
Circle(6)={3,1,7};
Circle(7)={7,1,5};
Circle(8)={5,1,6};

Line Loop(1)={-1,-6,-7,-4};
Line Loop(2)={1,-5,-8,4};
Line Loop(3)={-2,-3,7,6};
Line Loop(4)={2,3,8,5};

Ruled Surface(1)={1} In Sphere{1};
Ruled Surface(2)={2} In Sphere{1};
Ruled Surface(3)={3} In Sphere{1};
Ruled Surface(4)={4} In Sphere{1};
"""

    fn = tempname()
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname()
    run(`gmsh $fn -2 -format msh -o $fno`)
    fdo = open(fno,"r")
    m = read_gmsh_mesh(fdo)
    close(fdo)
    rm(fno)
    rm(fn)

    return m

end
