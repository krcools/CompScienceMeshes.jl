using Documenter, CompScienceMeshes

makedocs(;
    modules=[CompScienceMeshes],
    authors="Kristof Cools and contributors",
    sitename="CompScienceMeshes.jl Docs",
    checkdocs=:none,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://krcools.github.io/CompScienceMeshes.jl",
        edit_link="master",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=true,
    ))

# makedocs(
#     clean=false,
#     sitename = "CompScienceMeshes.jl Manual")
# # deploydocs(
# #     julia = "nightly",
# #     deps = Deps.pip("mkdocs", "python-markdown-math"),
# #     repo = "github.com/krcools/CompScienceMeshes.jl.git",)
# deploydocs(
#     repo = "github.com/krcools/CompScienceMeshes.jl.git",
# )

deploydocs(;
    repo="github.com/krcools/CompScienceMeshes.jl.git",
    target="build",
    push_preview=true,
    forcepush=true,
    #devbranch = "feature/docs",
)