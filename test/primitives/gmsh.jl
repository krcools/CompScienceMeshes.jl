@testitem "meshgeo [p=1,dim=2]" begin
    
    fn = joinpath(pkgdir(CompScienceMeshes), "test", "assets", "cuboid.geo")

    m = meshgeo(fn; element=:triangle)
    @test m isa CompScienceMeshes.Mesh
    @test dimension(m) == 2
end


@testitem "meshgeo [p=2,dim=2]" begin
    
    fn = joinpath(pkgdir(CompScienceMeshes), "test", "assets", "cuboid.geo")

    m = meshgeo(fn; element=:triangle, order=2)
    @test m isa CompScienceMeshes.CurvilinearMesh{3,3,6,Float64,2}
    @test dimension(m) == 2
end