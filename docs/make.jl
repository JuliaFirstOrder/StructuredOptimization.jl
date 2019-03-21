using Documenter, StructuredOptimization, LinearAlgebra, DSP, FFTW, AbstractOperators

makedocs(
  modules = [StructuredOptimization],
  format = Documenter.HTML(),
  sitename = "StructuredOptimization",
  authors = "Niccolò Antonello and Lorenzo Stella",
  pages = [
  "Home"                  => "index.md",
  "Quick Tutorial Guide"  => "tutorial.md",
  "Expressions"           => "expressions.md",
  "Functions"             => "functions.md",
  "Solvers"               => "solvers.md",
  "Demos"                 => "demos.md",
  ],
)

deploydocs(
  repo   = "github.com/kul-forbes/StructuredOptimization.jl.git",
  target = "build",
)
