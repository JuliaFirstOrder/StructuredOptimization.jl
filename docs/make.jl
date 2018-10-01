using Documenter, StructuredOptimization, LinearAlgebra, DSP, AbstractFFTs, FFTW, AbstractOperators

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
  "Demos"                 => "demos.md",
  ],
)

deploydocs(
  repo   = "github.com/kul-forbes/StructuredOptimization.jl.git",
  julia  = "1.0",
  osname = "linux",
  target = "build",
  deps = nothing,
  make = nothing,
)
