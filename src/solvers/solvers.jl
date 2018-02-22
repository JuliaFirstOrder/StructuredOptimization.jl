
include("terms_extract.jl")
include("terms_properties.jl")
include("terms_splitting.jl")

abstract type Solver end

abstract type ForwardBackwardSolver <: Solver end


################################################################################
export PG, FPG

struct PG <: ForwardBackwardSolver
    kwargs::Array
    function PG(; kwargs...)
        new(kwargs)
    end
end

function FPG(; kwargs...)
    return PG(; kwargs..., fast=true)
end

function build_iterator(x, solver::PG; kwargs...)
    x, ProximalAlgorithms.FBSIterator(~x; solver.kwargs..., kwargs...)
end

################################################################################
export ZeroFPR

struct ZeroFPR <: ForwardBackwardSolver
    kwargs::Array
    function ZeroFPR(; kwargs...)
        new(kwargs)
    end
end

function build_iterator(x, solver::ZeroFPR; kwargs...)
    x, ProximalAlgorithms.ZeroFPRIterator(~x; solver.kwargs..., kwargs...)
end

################################################################################
export PANOC

struct PANOC <: ForwardBackwardSolver
    kwargs::Array
    function PANOC(; kwargs...)
        new(kwargs)
    end
end

function build_iterator(x, solver::PANOC; kwargs...)
    x, ProximalAlgorithms.PANOCIterator(~x; solver.kwargs..., kwargs...)
end

default_solver = PANOC

################################################################################
export build

function build(terms::Tuple, solver::ForwardBackwardSolver)
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
			if is_linear(smooth)
				append!(kwargs, [(:As, As)])
			else
				fs = PrecomposeNonlinear(fs, As)
			end
			append!(kwargs, [(:fs, fs)])
		end
		return build_iterator(x, solver; kwargs...)
	end
	error("Sorry, I cannot solve this problem")
end

################################################################################
export solve!

function solve!(x_and_iter::Tuple{Tuple{Vararg{Variable}}, ProximalAlgorithms.ProximalAlgorithm})
    x, iterator = x_and_iter
    it, x_star = ProximalAlgorithms.run!(iterator)
    blockset!(~x, x_star)
    return iterator, it
end


export solve
function solve(terms::Tuple, solver::ForwardBackwardSolver)
    built_slv = build(terms, solver)
    solve!(built_slv)
end
