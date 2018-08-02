
#include("SparseDeconvolution.jl")
#results = SparseDeconvolution.run_demo()
#SparseDeconvolution.show_results(results...)
#
#results = SparseDeconvolution.run_demo_JuMP()
#SparseDeconvolution.show_results(results...)

#include("LineSpectraEstimation.jl")
#results = LineSpectraEstimation.run_demo()
###results = LineSpectraEstimation.run_demo_Convex()
#LineSpectraEstimation.show_results(results...)

#include("MatrixDecomposition.jl")
#results = MatrixDecomposition.run_demo()
#MatrixDecomposition.show_results(results...)
		
include("DNN.jl")
results = DNN.run_demo()
DNN.show_results(results...)

#include("TotalVariation.jl")
#results = TotalVariation.run_demo()
#TotalVariation.show_results(results...)

#include("AudioDeclipping.jl")
#results = AudioDeclipping.run_demo()
#AudioDeclipping.show_results(results...)
#AudioDeclipping.save_wav(results...)
