function solve!{N}(terms::NTuple{N,Term}, solver::PG)
	xAll = extract_variables(terms) 
	# Separate smooth and nonsmooth
	smooth, nonsmooth = split_smooth(terms)
	if is_proximable(xAll,nonsmooth)
		solver.verbose == true && println("--------------------Solving Primal---------------")
		f = extract_functions(smooth) 
		L = extract_operators(xAll,smooth) 
		g = extract_proximable(xAll,nonsmooth) 
		apply!(solver, ~xAll, f, L, g)
		return solver
	end
	#strongly = [t for t in terms if is_strongly_convex(t) == true]
	#nonstrongly = [t for t in terms if is_strongly_convex(t) == false]
	#if false # TODO: here, a condition for "easily conjugable" should go
	#	# Solving the DUAL
        #		solver.verbose == true && println("-------------------- Solving Dual ---------------")
	#	return solver
	#end
	error("Sorry, I cannot solve this problem")
end
