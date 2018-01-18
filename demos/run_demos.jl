
include("SparseDeconvolution.jl")
results = SparseDeconvolution.run_demo()
results = SparseDeconvolution.run_demo_cvx()
SparseDeconvolution.show_results(results...)

#include("LineSpectraEstimation.jl")
#results = LineSpectraEstimation.run_demo()
#results = LineSpectraEstimation.run_demo_cvx()
#LineSpectraEstimation.show_results(results...)
		
