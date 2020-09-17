using PCquery
using Documenter

makedocs(;
    modules=[PCquery],
    authors="Matt Karikomi <mattkarikomi@gmail.com> and contributors",
    repo="https://github.com/mkarikom/PCquery.jl/blob/{commit}{path}#L{line}",
    sitename="PCquery.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mkarikom.github.io/PCquery.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mkarikom/PCquery.jl",
)
