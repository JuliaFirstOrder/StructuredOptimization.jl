using Pkg
Pkg.activate(".")
Pkg.instantiate()

using IJulia
notebook(dir=pwd())
