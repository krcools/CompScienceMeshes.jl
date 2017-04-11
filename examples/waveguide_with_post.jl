using BoundaryElements
using CompScienceMeshes
using LinearForms

export meshwaveguidepost

"""
    meshwaveguidepost(width, height, length, post_radius, delta)

Creates a mesh for a waveguide of width (along the x-axis) `width`, height (along
    the y-axis) `height` and length (along the z-axis) `length` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

A cylindrical post of radius `post_radius` is positioned vertically (i.e.
    parallel to y-axis) midway along the waveguide.

The target edge size is `delta`.
"""

function meshwaveguidepost(width, height, length, post_radius, delta)
    s =
"""

lc = $delta;

//--------------WAVEGUIDE POINTS--------------

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.0, $height, 0.0, lc};
Point(3) = {$width, $height, 0.0, lc};
Point(4) = {$width, 0.0, 0.0, lc};

Point(5) = {0.0, 0.0, $length, lc};
Point(6) = {0.0, $height, $length, lc};
Point(7) = {$width, $height, $length, lc};
Point(8) = {$width, 0.0, $length, lc};


//--------------POST POINTS-------------------


Point(9) = {($width)/2,0,($length)/2,lc};

Point(10) = {($width)/2,0,(($length)/2)-(($post_radius)/2),lc};
Point(11) = {(($width)/2)-(($post_radius)/2),0,($length)/2,lc};
Point(12) = {($width)/2,0,(($length)/2)+(($post_radius)/2),lc};
Point(13) = {(($width)/2)+(($post_radius)/2),0,($length)/2,lc};

Point(14) = {($width)/2,$height,($length)/2,lc};

Point(15) = {($width)/2,$height,(($length)/2)-(($post_radius)/2),lc};
Point(16) = {(($width)/2)-(($post_radius)/2),$height,($length)/2,lc};
Point(17) = {($width)/2,$height,(($length)/2)+(($post_radius)/2),lc};
Point(18) = {(($width)/2)+(($post_radius)/2),$height,($length)/2,lc};

//--------------WAVEGUIDE-POST CONTINUOUS SURFACE DEFINITION---------

//Part 1 - waveguide boundaries

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,5};

Line(9)={1,5};
Line(10)={2,6};
Line(11)={3,7};
Line(12)={4,8};

//Part 2 - discontinuity boundaries

Circle(13) = {10,9,11};
Circle(14) = {11,9,12};
Circle(15) = {12,9,13};
Circle(16) = {13,9,10};

Circle(17) = {15,14,16};
Circle(18) = {16,14,17};
Circle(19) = {17,14,18};
Circle(20) = {18,14,15};

Line(21) = {10,15};
Line(22) = {11,16};
Line(23) = {12,17};
Line(24) = {13,18};

//Part 3 - creating skeleton

Line Loop(1)={-1,9,5,-10};

Line Loop(2)={-2,10,6,-11};
Line Loop(3)={17,18,19,20};

Line Loop(4)={-3,11,7,-12};

Line Loop(5)={-4,12,8,-9};
Line Loop(6)={-16,-15,-14,-13};

Line Loop(7) = {13,22,-17,-21};
Line Loop(8) = {14,23,-18,-22};
Line Loop(9) = {15,24,-19,-23};
Line Loop(10) = {16,21,-20,-24};

//Part 4 - creating surfaces

Plane Surface(1)={1};
Plane Surface(2)={3,2};
Plane Surface(3)={4};
Plane Surface(4)={6,5};

Ruled Surface(5) = {7};
Ruled Surface(6) = {8};
Ruled Surface(7) = {9};
Ruled Surface(8) = {10};

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
    # include(Pkg.dir("CompScienceMeshes","src","gmsh.jl"))
    m = CompScienceMeshes.read_gmsh_mesh(fdo)
    # close(fdo)
    # rm(fno)
    # rm(fn)

    return m

end

# Γ = meshwaveguidepost(22,22,50,3.81,0.8)
# RT = raviartthomas(Γ)
#
# include(Pkg.dir("CompScienceMeshes","examples","plotlyjs_patches.jl"))
# A = rand(numcells(Γ))
# p = patch(Γ, A)
# PlotlyJS.plot([p])
