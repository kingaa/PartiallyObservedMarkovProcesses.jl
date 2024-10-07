push!(LOAD_PATH,"../src/")

using POMP
using Documenter

makedocs(
    sitename = "POMP.jl",
    modules  = [POMP],
    repo = Remotes.GitHub("kingaa","POMP.jl"),
    pages=[
        "Home" => "index.md"
    ]
)

deploydocs(
    ;
    repo="github.com/kingaa/POMP.jl"
)
