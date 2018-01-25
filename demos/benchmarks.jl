using BenchmarkTools
BLAS.set_num_threads(4)

#verbose, samples, seconds
v,        smp,     sec      = 0, 5, 100
suite = BenchmarkGroup()

suite["SparseDeconvolution"] = BenchmarkGroup()
include("SparseDeconvolution.jl")
suite["SparseDeconvolution"]["FullMatrix"] = SparseDeconvolution.benchmark(verb = v, samples = smp, seconds = sec)
suite["SparseDeconvolution"]["MatrixFree"] = SparseDeconvolution.benchmarkMatrixFree(verb = v, samples = smp, seconds = sec)

suite["LineSpectraEstimation"] = BenchmarkGroup()
include("LineSpectraEstimation.jl")
suite["LineSpectraEstimation"]["FullMatrix"] = LineSpectraEstimation.benchmark(verb = v, samples = smp, seconds = sec)
suite["LineSpectraEstimation"]["MatrixFree"] = LineSpectraEstimation.benchmarkMatrixFree(verb = v, samples = smp, seconds = sec)

suite["MatrixDecomposition"] = BenchmarkGroup()
include("MatrixDecomposition.jl")
suite["MatrixDecomposition"] = MatrixDecomposition.benchmark(verb = v, samples = smp, seconds = sec)

suite["TotalVariation"] = BenchmarkGroup()
include("TotalVariation.jl")
suite["TotalVariation"] = TotalVariation.benchmark(verb = v, samples = smp, seconds = sec)

suite["DNN"] = BenchmarkGroup()
include("DNN.jl")
suite["DNN"] = DNN.benchmark(verb = v, samples = smp, seconds = sec)

suite["TotalVariation"] = BenchmarkGroup()
include("TotalVariation.jl")
suite["TotalVariation"] = TotalVariation.benchmark(verb = v, samples = smp, seconds = sec)

suite["AudioDeclipping"] = BenchmarkGroup()
include("AudioDeclipping.jl")
suite["AudioDeclipping"] = AudioDeclipping.benchmark(verb = v, samples = 1, seconds = sec)

println("\n")
showall(median(suite))
println("\n")
#showall(memory(suite))
