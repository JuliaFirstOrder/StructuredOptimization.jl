export solve!

abstract type Solver end

include("terms_extract.jl")
include("terms_properties.jl")
include("terms_splitting.jl")

include("fbsolvers/fbsolvers.jl")
include("fbsolvers/pg.jl")
include("fbsolvers/zerofpr.jl")

default_slv = ZeroFPR
