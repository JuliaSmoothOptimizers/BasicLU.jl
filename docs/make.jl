using Documenter, BasicLU

makedocs(
  modules = [BasicLU],
  doctest = true,
  linkcheck = true,
  strict = true,
  format = Documenter.HTML(assets = ["assets/style.css"],
                           ansicolor = true,
                           prettyurls = get(ENV, "CI", nothing) == "true"),
  sitename = "BasicLU.jl",
  pages = ["Home" => "index.md"],
)

deploydocs(
  repo = "github.com/JuliaSmoothOptimizers/BasicLU.jl.git",
  push_preview = true,
  devbranch = "main",
)
