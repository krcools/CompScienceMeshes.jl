@testitem "load_gmsh_mesh: physical group spanning multiple entities" begin

    fn = joinpath(pkgdir(CompScienceMeshes), "test", "assets", "cuboid_triangle.msh")

    # "Side" consists of four geometrical surfaces, "Top" and "Bottom" of one each.
    top    = load_gmsh_mesh(fn, physical="Top")
    bottom = load_gmsh_mesh(fn, physical="Bottom")
    side   = load_gmsh_mesh(fn, physical="Side")
    whole  = load_gmsh_mesh(fn)

    # All elements of all four surfaces must be loaded, not just those of the first.
    @test numcells(side) == 4 * numcells(top)
    @test numcells(top) + numcells(bottom) + numcells(side) == numcells(whole)

    # Agrees with the legacy reader, which filters elements on the physical tag.
    side_old = CompScienceMeshes.read_gmsh_mesh(fn, physical="Side")
    @test numcells(side) == numcells(side_old)
    sorted_cells(m) = sort([sort(c.indices) for c in cells(m)])
    @test sorted_cells(side) == sorted_cells(side_old)

    # The four sides of a cube form a tube: its boundary consists of two square loops.
    @test CompScienceMeshes.isoriented(side)
    @test numcells(boundary(side)) == numcells(boundary(top)) + numcells(boundary(bottom))
end
