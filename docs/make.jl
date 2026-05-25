using Documenter
using GMRF

makedocs(
    sitename = "GMRF.jl",
    modules = [GMRF],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/jojal5/GMRF.jl.git",
    devbranch = "master",
)
