export build

"""
`build(terms::Tuple, solver_opt::ForwardBackwardSolver)`

Takes as input a tuple containing the terms defining the problem and the solver options.

Returns a tuple containing the optimization variables and the built solver.

# Example

```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem( ls(A*x - b ) , norm(x) <= 1 );

julia> build(p, PG());

```

"""
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
      if is_linear(smooth)
        fs = extract_functions(smooth)
        As = extract_operators(x, smooth)
        append!(kwargs, [(:As, As)])
      else
        fs = extract_functions_nodisp(smooth)
        As = extract_affines(x, smooth)
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

"""
`solve!( x_solver )`

Takes as input a tuple containing the optimization variables and the built solver.

Solves the problem returning a tuple containing the iterations taken and the build solver.

# Example

```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> p = problem( ls(A*x - b ) , norm(x) <= 1 );

julia> x_solver = build(p, PG(verbose = 0));

julia> solve!(x_solver);

```

"""
function solve!(x_and_iter::Tuple{Tuple{Vararg{Variable}}, ProximalAlgorithms.ProximalAlgorithm})
  x, iterator = x_and_iter
  it, x_star = ProximalAlgorithms.run!(iterator)
  ~x .= x_star
  return it, iterator 
end


export solve

"""
`solve(terms::Tuple, solver_opt::ForwardBackwardSolver)`

Takes as input a tuple containing the terms defining the problem and the solver options.

Solves the problem returning a tuple containing the iterations taken and the build solver.

# Example

```julia

julia> x = Variable(4)
Variable(Float64, (4,))

julia> A, b = randn(10,4), randn(10);

julia> solve(p,PG());
it |      gamma |        fpr |
------|------------|------------|
1 | 7.6375e-02 | 1.8690e+00 |
12 | 7.6375e-02 | 9.7599e-05 |

```

"""
function solve(terms::Tuple, solver::ForwardBackwardSolver)
  built_slv = build(terms, solver)
  return solve!(built_slv)
end
