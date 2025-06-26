using JSOTemplate
using Documenter

DocMeta.setdocmeta!(BilevelBenchmark, :DocTestSetup, :(using BilevelBenchmark); recursive = true)

makedocs(;
  modules = [BilevelBenchmark],
  doctest = true,
  linkcheck = false,
  strict = false,
  authors = "Valentin Dijon <vanadjy@gmail.com> and contributors",
  repo = "https://github.com/JuliaSmoothOptimizers/BilevelBenchmark.jl/blob/{commit}{path}#{line}",
  sitename = "BilevelBenchmark.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://github.com/Vanadjy/BilevelBenchmark.jl",
    assets = ["assets/style.css"],
  ),
  pages = ["Home" => "index.md", "Reference" => "reference.md"],
)

deploydocs(;
  repo = "github.com/Vanadjy/BilevelBenchmark.jl",
  push_preview = true,
  devbranch = "main",
)
