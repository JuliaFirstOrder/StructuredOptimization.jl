using Documenter, StructuredOptimization

makedocs(
  modules = [StructuredOptimization],
  format = :html,
  sitename = "StructuredOptimization",
  authors = "Niccolò Antonello and Lorenzo Stella",
  pages = Any[
  "Home"            => "index.md",
  "Variables"       => "variables.md",
  "Expressions"     => "expressions.md",
  ],
)

deploydocs(
  repo   = "github.com/kul-forbes/StructuredOptimization.jl.git",
  julia  = "0.6",
  osname = "linux",
  target = "build",
  deps = nothing,
  make = nothing,
)
