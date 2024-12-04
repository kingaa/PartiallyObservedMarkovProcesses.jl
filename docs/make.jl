using Documenter, POMP, RCall

makedocs(
    sitename = "POMP.jl",
    modules  = [POMP],
    repo = Remotes.GitHub("kingaa","POMP.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        size_threshold = 600 * 2^10,
        size_threshold_warn = 500 * 2^10, # 600 KiB
        canonical="https://github.com/kingaa/POMP.jl/",
        edit_link="master",
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(
    ;
    repo="github.com/kingaa/POMP.jl"
)
