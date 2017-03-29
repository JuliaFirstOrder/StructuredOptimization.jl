abstract ExtendedRealValuedFunction
abstract SmoothFunction      <: ExtendedRealValuedFunction
abstract NonSmoothFunction   <: ExtendedRealValuedFunction

include("functions/CostFunction.jl")
include("functions/absorb_merge.jl")
include("functions/LeastSquares.jl")
include("functions/Norm.jl")
include("functions/Box.jl")
include("functions/HingeLoss.jl")

gradient(f::ExtendedRealValuedFunction) = error("gradient not implemented for $f")

function Base.show(io::IO, cf::CostFunction)
	if isempty(cf)
		println(io, "Empty Cost Function") 
	else
		description = fun_name(cf.f[1],1)
		operator    = "\n A1 = "*fun_name(RegLS.operator(cf.A[1]))
		parameter   = fun_par(cf.f[1],1)
		for i = 2:length(cf.f)
			description = description*"+ "fun_name(cf.f[i],i)
			operator    = operator*",\n A$i = "*fun_name(RegLS.operator(cf.A[i]))
			parameter = parameter*", "fun_par(cf.f[i],i)
		end
		
		println(io, "description : ", description) 
		println(io, "operators   : ", operator   )
		println(io, "parameters  : ", parameter  )
	end
end





