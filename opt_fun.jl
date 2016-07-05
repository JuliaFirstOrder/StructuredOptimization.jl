function prox_l1(x,λ)

	#prox operator of g(x) =  || x ||_1
	return max(0, x - λ) - max(0, -x - λ)

end


function fista(A,b,λ,x0,MaxIt,TOL,prox_fun)

	##############
	# Inputs:
	# A Matrix ∈ R^{m ⋅ n} 
	# b vector ∈ R^{m}
	# x0 initial guess, vector ∈ R^{n}
	# MaxIt maximum number of iterations
	# TOL tolerance
	# prox_fun function of a proximal operator 
	##############
	# Outputs:
	# x0 solution
	# cf vector containing cost function value
	#    as a function of iteration k
	##############

	α = 1    # step-size
	β = 0.5  # parameter used to decrease step-size

	#vectors used to accelerate the method
	z = zeros(x0)
	x_prev = copy(x0)
	y = zeros(x0)

	cf = zeros(MaxIt) #initialize cost function vector 

	for k = 1:MaxIt
	
		
                cf[k] = 0.5*norm(A*x0-b)
		y = x0+(k/(k+3))*(x0-x_prev)
                grad_y = A'*(A*y-b)

		while(true)
	
			z = prox_fun(y-α*grad_y, α*λ)
		
			if(( 0.5*norm(A*x-b) #linesearch
		            -0.5*norm(A*y-b)-grad_y'*(z-y)-(1/(2*α))*norm(z-y)^2)[1]<=0  )
				break
			end
			α = β*α #decrease step-size
		end


		residual = norm(z-x0) #compute residual 
		if(residual<=TOL)
			print_stuff(k,MaxIt,residual)
			println()
			println("acceptable solution found")
			cf = cf[1:k] #cut vector
			break
		end

		print_stuff(k,MaxIt,residual)
		x_prev = copy(x0)
		x0 = copy(z)
		
	end
	println()

	return x0,cf


end


function print_stuff(k,MaxIt,residual)

	print("\u1b[1G")   # go to first column
	print("it: ",k," / ",MaxIt,"  Res: ", residual)
	print("\u1b[K")    # clear the rest of the line
end
