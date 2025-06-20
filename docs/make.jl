using IdentityByDescentDispersal
using Documenter

DocMeta.setdocmeta!(
    IdentityByDescentDispersal,
    :DocTestSetup,
    :(using IdentityByDescentDispersal);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [IdentityByDescentDispersal],
    authors = "Curro Campuzano Jiménez campuzanocurro@gmail.com",
    repo = "https://github.com/currocam/IdentityByDescentDispersal.jl/blob/{commit}{path}#{line}",
    sitename = "IdentityByDescentDispersal.jl",
    format = Documenter.HTML(;
        canonical = "https://currocam.github.io/IdentityByDescentDispersal.jl",
    ),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/currocam/IdentityByDescentDispersal.jl")
