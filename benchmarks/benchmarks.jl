using BenchmarkTools
BLAS.set_num_threads(4)

verb, samples, seconds = 0, 5, 60^2
suite = BenchmarkGroup()

suite["SparseDeconvolution"] = BenchmarkGroup()
include("SparseDeconvolution.jl")
suite["SparseDeconvolution"]["FullMatrix"] = SparseDeconvolution.benchmark(verb = verb, samples = samples, seconds = seconds)
suite["SparseDeconvolution"]["MatrixFree"] = SparseDeconvolution.benchmarkMatrixFree(verb = verb, samples = samples, seconds = seconds)

suite["LineSpectraEstimation"] = BenchmarkGroup()
include("LineSpectraEstimation.jl")
suite["LineSpectraEstimation"]["FullMatrix"] = LineSpectraEstimation.benchmark(verb = verb, samples = samples, seconds = seconds)
suite["LineSpectraEstimation"]["MatrixFree"] = LineSpectraEstimation.benchmarkMatrixFree(verb = verb, samples = samples, seconds = seconds)

suite["MatrixDecomposition"] = BenchmarkGroup()
include("MatrixDecomposition.jl")
suite["MatrixDecomposition"] = MatrixDecomposition.benchmark(verb = verb, samples = samples, seconds = seconds)

suite["DNN"] = BenchmarkGroup()
include("DNN.jl")
suite["DNN"] = DNN.benchmark(verb = verb, samples = samples, seconds = seconds)

suite["TotalVariation"] = BenchmarkGroup()
include("TotalVariation.jl")
suite["TotalVariation"] = TotalVariation.benchmark(verb = verb, samples = samples, seconds = seconds)

suite["AudioDeclipping"] = BenchmarkGroup()
include("AudioDeclipping.jl")
suite["AudioDeclipping"] = AudioDeclipping.benchmark(verb = verb, samples = 1, seconds = seconds)

println("\n")
showall(median(suite))
println("\n")
#showall(memory(suite))
