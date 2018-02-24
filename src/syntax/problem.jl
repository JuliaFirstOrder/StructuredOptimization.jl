export problem

"""
`problems(terms...)`

Constructs a problem.

# Example

```julia

julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem(ls(A*x-b), norm(x) <= 1)

```

"""
function problem(terms::Vararg)
	cf = ()
	for i = 1:length(terms)
		cf = (cf...,terms[i]...)
	end
	return cf
end
