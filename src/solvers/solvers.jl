export solve!

include("terms_extract.jl")
include("terms_properties.jl")
include("terms_splitting.jl")

abstract type Solver end

abstract type FBSolver <: Solver end

################################################################################

struct PG <: FBSolver
    kwargs::Array
    function PG(; kwargs...)
        new(kwargs)
    end
end

function FPG(; kwargs...)
    return PG(; kwargs..., fast=true)
end

function apply!(solver::PG, x; kwargs...)
    (it, xsol, solver) = ProximalAlgorithms.FBS(x; solver.kwargs..., kwargs...)
    blockcopy!(x, xsol)
    return solver
end

################################################################################

struct ZeroFPR <: FBSolver
    kwargs::Array
    function ZeroFPR(; kwargs...)
        new(kwargs)
    end
end

function apply!(solver::ZeroFPR, x; kwargs...)
    (it, xsol, solver) = ProximalAlgorithms.ZeroFPR(x; solver.kwargs..., kwargs...)
    blockcopy!(x, xsol)
    return solver
end

################################################################################

function solve!(terms::Tuple, solver::FBSolver)
	x = extract_variables(terms)
	# Separate smooth and nonsmooth
	smooth, nonsmooth = split_smooth(terms)
	# Separate quadratic and nonquadratic
	quadratic, smooth = split_quadratic(smooth)
    kwargs = Array{Any, 1}()
	if is_proximable(nonsmooth)
        g = extract_proximable(x, nonsmooth)
        append!(kwargs, [(:g, g)])
		if !isempty(quadratic)
			fq = extract_functions(quadratic)
			Aq = extract_operators(x, quadratic)
            append!(kwargs, [(:fq, fq)])
            append!(kwargs, [(:Aq, Aq)])
        end
		if !isempty(smooth)
			fs = extract_functions(smooth)
			As = extract_operators(x, smooth)
            append!(kwargs, [(:fs, fs)])
            append!(kwargs, [(:As, As)])
		end
		return apply!(solver, ~x; kwargs...)
	end
	error("Sorry, I cannot solve this problem")
end

################################################################################

default_solver = ZeroFPR
