using Documenter
using ParticleTracer

DocMeta.setdocmeta!(ParticleTracer, :DocTestSetup, :(using ParticleTracer); recursive=true)


push!(LOAD_PATH,"../src/")
makedocs(sitename="ParticleTracer.jl",
         doctest = false,
         pages = [
            "Home" => "index.md",
            "Examples" => "examples.md",
            "API reference" => "api.md",
            "Notes" => "notes.md",
            "Contact Developer" => "Contact.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/junyixu/ParticleTracer.jl.git",
    devbranch = "TestPtcs/doc"
)
