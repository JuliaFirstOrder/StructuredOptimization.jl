using BenchmarkTools

#verbose, samples, seconds
v,        smp,     sec      = 0, 10, 100
suite = BenchmarkGroup()

suite["SparseDeconvolution"] = BenchmarkGroup()
include("SparseDeconvolution.jl")
suite["SparseDeconvolution"]["FullMatrix"] = SparseDeconvolution.benchmark(verb = v, samples = smp, seconds = sec)
suite["SparseDeconvolution"]["MatrixFree"] = SparseDeconvolution.benchmarkMatrixFree(verb = v, samples = smp, seconds = sec)

suite["LineSpectraEstimation"] = BenchmarkGroup()
include("LineSpectraEstimation.jl")
suite["LineSpectraEstimation"]["FullMatrix"] = LineSpectraEstimation.benchmark(verb = v, samples = smp, seconds = sec)
suite["LineSpectraEstimation"]["MatrixFree"] = LineSpectraEstimation.benchmarkMatrixFree(verb = v, samples = smp, seconds = sec)

println("\n")
showall(median(suite))
println("\n")
#showall(memory(suite))
