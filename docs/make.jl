using IdentityByDescentDispersal
using Documenter
using Literate

# Set up Documenter doctest environment
DocMeta.setdocmeta!(
    IdentityByDescentDispersal,
    :DocTestSetup,
    :(using IdentityByDescentDispersal);
    recursive = true,
)

# Generate tutorial markdown from Literate source
Literate.markdown("docs/src/tutorial.jl", "docs/src/"; name = "tutorial", documenter = true)
Literate.markdown("docs/src/overview.jl", "docs/src/"; name = "overview", documenter = true)

# Optional: Rename pages for display
const page_rename = Dict("developer.md" => "Developer docs")

# Collect all markdown files except index.md and tutorial.md
const numbered_pages = sort([
    file for file in readdir("docs/src") if
    endswith(file, ".md") && file ∉ ("index.md", "tutorial.md", "overview.md")
])

makedocs(;
    modules = [IdentityByDescentDispersal],
    authors = "Curro Campuzano Jiménez <campuzanocurro@gmail.com>",
    repo = "https://github.com/currocam/IdentityByDescentDispersal.jl/blob/{commit}{path}#{line}",
    sitename = "IdentityByDescentDispersal.jl",
    format = Documenter.HTML(;
        canonical = "https://currocam.github.io/IdentityByDescentDispersal.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Theory overview" => "overview.md",
        "Developer docs" => numbered_pages,
    ],
)

deploydocs(; repo = "github.com/currocam/IdentityByDescentDispersal.jl")
