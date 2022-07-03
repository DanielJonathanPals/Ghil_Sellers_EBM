using Ghil_Sellers_EBM
using Documenter

DocMeta.setdocmeta!(Ghil_Sellers_EBM, :DocTestSetup, :(using Ghil_Sellers_EBM); recursive=true)

makedocs(;
    modules=[Ghil_Sellers_EBM],
    authors="Daniel Pals <Daniel.Pals@tum.de>",
    repo="https://github.com/DanielJonathanPals/Ghil_Sellers_EBM.jl/blob/{commit}{path}#{line}",
    sitename="Ghil_Sellers_EBM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
