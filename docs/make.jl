using Documenter, CompScienceMeshes

makedocs(clean=false)
deploydocs(
    julia = "nightly",
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/krcools/CompScienceMeshes.jl.git",
)
