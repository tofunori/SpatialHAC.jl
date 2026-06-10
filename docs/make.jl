using Documenter, SpatialHAC

makedocs(;
    modules = [SpatialHAC],
    sitename = "SpatialHAC.jl",
    authors = "Thierry Laurent St-Pierre",
    format = Documenter.HTML(;
        canonical = "https://tofunori.github.io/SpatialHAC.jl",
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "The estimator" => "estimator.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/tofunori/SpatialHAC.jl",
    devbranch = "main",
    push_preview = false,
)
