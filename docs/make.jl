using Documenter, StructuredOptimization

makedocs(
  modules = [StructuredOptimization],
  format = :html,
  sitename = "StructuredOptimization",
  authors = "Niccolò Antonello and Lorenzo Stella",
  pages = Any[
  "Home"                  => "index.md",
  "Quick Tutorial Guide"  => "tutorial.md",
  "Expressions"           => "expressions.md",
  "Functions"             => "functions.md",
  "Solvers"               => "solvers.md",
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
