using GLVisualize, GeometryTypes, Colors, GLAbstraction
using CompScienceMeshes

const cm = vec(mapslices(
    r->Colors.RGBA{Float32}(r[1],r[2],r[3], 1.0),
    readcsv(Pkg.dir("CompScienceMeshes","examples","cm.csv")),
    [2]))

function patch(mesh::Mesh, amp = nothing)

    if amp == nothing
        amp = Float32[v[end] for v in mesh.vertices]
    else
        vtc, nc = vertextocellmap(mesh)
        amp = Float32[mean(amp[vtc[i,1:nc[i]]]) for i in 1:numvertices(mesh)]
    end

    @show length(amp)

    m, M = extrema(amp)

    #cm = map(x -> RGBA{Float32}(x, 1.0), colormap("RdBu"))
    nc = length(cm)
    at = Float32[(nc-1)*(a-m)/(M-m) for a in amp]

    @assert length(at) == numvertices(mesh)

    vertices = Point3f0[v for v in mesh.vertices]
    faces = GLTriangle[f .- 1 for f in mesh.faces]

    #colors = RGBA{U8}[RGBA{U8}(rand(), rand(), rand(), 1.) for i=1:5]
    #attribute_id = rand(0f0:4f0, length(faces))

    sphere_mesh = GLNormalAttributeMesh(
        vertices=vertices, faces=faces,
        attributes=cm, attribute_id=at)

end

## Example
# m = meshsphere(1.0, 0.2)
# sphere_mesh = patch(m)
#
error("stop")

##
fcr, geo = facecurrents(u, RT)
mp = patch(geo,real.(norm.(fcr)))


##
window = glscreen()
moveright = translationmatrix(Vec3f0(0,0,0))
_view(visualize(mp, model=moveright), window)
renderloop(window)
