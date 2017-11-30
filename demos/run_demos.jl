

#include("SparseDeconvolution.jl")
#results = SparseDeconvolution.run_demo()
#SparseDeconvolution.show_results(results...)

#results = SparseDeconvolution.run_demo_Convex()
#SparseDeconvolution.show_results(results...)

#results = SparseDeconvolutionConvex.run_demo();
#SparseDeconvolution.show_results(results...)

include("LineSpectraEstimation.jl")

results = LineSpectraEstimation.run_demo()
LineSpectraEstimation.show_results(results...)

results1 = LineSpectraEstimation.run_demo_Convex()
