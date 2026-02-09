using IdentityByDescentDispersal
using Documenter
using DocumenterCitations
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
Literate.markdown(
    "docs/src/inference.jl",
    "docs/src/";
    name = "inference",
    documenter = true,
)

# Optional: Rename pages for display
const page_rename = Dict("developer.md" => "Developer docs")

# Collect all markdown files except index.md and tutorial.md
const numbered_pages = sort([
    file for file in readdir("docs/src") if endswith(file, ".md") &&
    file ∉ ("index.md", "tutorial.md", "overview.md", "inference.md", "references.md")
])

# Set up bibliography
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules = [IdentityByDescentDispersal],
    authors = "Curro Campuzano Jiménez <campuzanocurro@gmail.com>",
    repo = "https://github.com/currocam/IdentityByDescentDispersal.jl/blob/{commit}{path}#{line}",
    sitename = "IdentityByDescentDispersal.jl",
    format = Documenter.HTML(;
        canonical = "https://currocam.github.io/IdentityByDescentDispersal.jl",
    ),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        "Theory overview" => "overview.md",
        "Basic usage" => "tutorial.md",
        "Inference and model evaluation" => "inference.md",
        "Developer docs" => numbered_pages,
        "References" => "references.md",
    ],
)

deploydocs(; repo = "github.com/currocam/IdentityByDescentDispersal.jl")
