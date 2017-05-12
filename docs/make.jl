using Documenter, CompScienceMeshes

makedocs(clean=false)
deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/krcools/CompScienceMeshes.jl.git",
)
