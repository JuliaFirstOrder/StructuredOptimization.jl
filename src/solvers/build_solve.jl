export build

"""
    parse_problem(terms::Tuple, solver::ForwardBackwardSolver)

Takes as input a tuple containing the terms defining the problem and the solver.

Returns a tuple containing the optimization variables and the problem terms
to be fed into the solver.

# Example

```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem( ls(A*x - b ) , norm(x) <= 1 );

julia> build(p, PG());
```
"""
function parse_problem(terms::Tuple, solver::T) where T <: ForwardBackwardSolver
  x = extract_variables(terms)
  # Separate smooth and nonsmooth
  smooth, nonsmooth = split_smooth(terms)
  if is_proximable(nonsmooth)
    g = extract_proximable(x, nonsmooth)
    kwargs = Dict{Symbol, Any}(:g => g)
    if !isempty(smooth)
      if is_linear(smooth)
        f = extract_functions(smooth)
        A = extract_operators(x, smooth)
        kwargs[:A] = A
      else  # ??
        f = extract_functions_nodisp(smooth)
        A = extract_affines(x, smooth)
        f = PrecomposeNonlinear(f, A)
      end
      kwargs[:f] = f
    end
    return (x, kwargs)
  end
  error("Sorry, I cannot parse this problem for solver of type $(T)")
end


export solve

"""
    solve(terms::Tuple, solver::ForwardBackwardSolver)

Takes as input a tuple containing the terms defining the problem and the solver options.

Solves the problem returning a tuple containing the iterations taken and the build solver.

# Example

```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem(ls(A*x - b ), norm(x) <= 1);

julia> solve(p, ProximalAlgorithms.ForwardBackward());

julia> ~x
4-element Array{Float64,1}:
 -0.6427139974173074
 -0.29043653211431103
 -0.6090539651510192
  0.36279278640995494
```
"""
function solve(terms::Tuple, solver::ForwardBackwardSolver)
    x, kwargs = parse_problem(terms, solver)
    x_star, it = solver(~x; kwargs...)
    ~x .= x_star
    return x, it
end
