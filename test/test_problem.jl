@printf("\nTesting extraction from Terms\n")

# testing extracting stuff from terms

m,n1 = 5,3
x1 = Variable(n1)
A = randn(m,n1)
# single term, single variable
cf = ls(A*x1)
xAll = RegLS.extract_variables(cf)
@test xAll[1] == x1
L = RegLS.extract_operators(xAll,cf)
@test typeof(L) <: MatrixOp
f = RegLS.extract_functions(cf)
@test typeof(f) <: SqrNormL2

# multiple terms, single variable
b1 = randn(n1)
cf = ls(A*x1) + 2.5*norm(x1+b1,1)
xAll = RegLS.extract_variables(cf)
@test xAll[1] == x1
V = RegLS.extract_operators(xAll,cf)
@test typeof(V) <: VCAT
@test typeof(V[1]) <: MatrixOp
@test typeof(V[2]) <: Eye
f = RegLS.extract_functions(cf)
@test typeof(f) <: SeparableSum
@test typeof(f.fs[1]) <: SqrNormL2
@test typeof(f.fs[2]) <: Postcompose
x = randn(n1)
@test norm(f.fs[2](x) - 2.5*norm(x+b1,1)) < 1e-12

# single term, multiple variables
x2 = Variable(m)
cf = ls(A*x1+x2+20)
xAll = RegLS.extract_variables(cf)
xAll = (x2,x1) # change the order on pourpose
H = RegLS.extract_operators(xAll,cf)
@test typeof(H) <: HCAT
@test typeof(H[1]) <: Eye
@test typeof(H[2]) <: MatrixOp
f = RegLS.extract_functions(cf)
@test typeof(f) <: PrecomposeDiagonal

## multiple terms, multiple variables
n1,n2,n3,n4,n5 = 3,3,4,4,7
A = randn(n5,n1)
x1,x2,x3,x4,x5 = Variable(randn(n1)),Variable(randn(n2)),Variable(randn(n3)),Variable(randn(n4)),Variable(randn(n5))

cf = ls(x1+x2)
xAll = RegLS.extract_variables(cf)
@test xAll == (x1,x2)

cf = ls(x1+x2)+ls(x1)
xAll = RegLS.extract_variables(cf)
@test xAll == (x1,x2)

cf = ls(x1+x2+3)+ls(x3+x4)+ls(x5)+ls(x5+A*x2)+ls(x1)+ls(x5)
xAll = RegLS.extract_variables(cf)
@test xAll == (x1,x2,x3,x4,x5)

V = RegLS.extract_operators(xAll,cf)

@test typeof(V[1][1]) <: Eye
@test typeof(V[1][2]) <: Eye
@test typeof(V[1][3]) <: Zeros
@test typeof(V[1][4]) <: Zeros
@test typeof(V[1][5]) <: Zeros

@test typeof(V[2][1]) <: Zeros
@test typeof(V[2][2]) <: Zeros
@test typeof(V[2][3]) <: Eye
@test typeof(V[2][4]) <: Eye
@test typeof(V[2][5]) <: Zeros

@test typeof(V[3][1]) <: Zeros
@test typeof(V[3][2]) <: Zeros
@test typeof(V[3][3]) <: Zeros
@test typeof(V[3][4]) <: Zeros
@test typeof(V[3][5]) <: Eye

@test typeof(V[4][1]) <: Zeros
@test typeof(V[4][2]) <: MatrixOp
@test typeof(V[4][3]) <: Zeros
@test typeof(V[4][4]) <: Zeros
@test typeof(V[4][5]) <: Eye

@test typeof(V[5][1]) <: Eye
@test typeof(V[5][2]) <: Zeros
@test typeof(V[5][3]) <: Zeros
@test typeof(V[5][4]) <: Zeros
@test typeof(V[5][5]) <: Zeros

@test typeof(V[6][1]) <: Zeros
@test typeof(V[6][2]) <: Zeros
@test typeof(V[6][3]) <: Zeros
@test typeof(V[6][4]) <: Zeros
@test typeof(V[6][5]) <: Eye

@printf("\nTesting splitting Terms\n")

x = Variable(5)
y = Variable(5)
cf = ls(x)+10*norm(x,2)+ls(x+y)

smooth, nonsmooth = RegLS.split_smooth(cf)
@test smooth[1] == cf[1]
@test smooth[2] == cf[3]
@test nonsmooth[1] == cf[2]

cf = ls(x)
smooth, nonsmooth = RegLS.split_smooth((cf,))
@test smooth == (cf,)
@test nonsmooth == ()

cf = norm(x,1)+norm(y,2)+norm(randn(5,5)*x+y,Inf)
xAll = RegLS.extract_variables(cf)
AAc, nonAAc = RegLS.split_AAc_diagonal(cf)
@test AAc[1] == cf[1]
@test AAc[2] == cf[2]
@test nonAAc[1] == cf[3]

cf = ls(sigmoid(x)) + ls(x)
quad, smooth = RegLS.split_quadratic(cf)
@test smooth[1] == cf[1]
@test quad[1] == cf[2]

@printf("\nTesting extracting Proximable functions\n")
# testing is_proximable
@test RegLS.is_proximable(AAc) == true
@test RegLS.is_proximable(nonAAc) == false

cf = norm(x[1:2],1)+norm(x[3:5])
xAll = RegLS.extract_variables(cf)

@test all(RegLS.is_AAc_diagonal.(cf)) == true
# @test RegLS.is_proximable(cf) == true

cf = norm(x[1:2],1)+norm(x[3:5])+norm(x,Inf)
xAll = RegLS.extract_variables(cf)

@test all(RegLS.is_AAc_diagonal.(cf)) == true
@test RegLS.is_proximable(cf) == false

# testing extract_proximable
# single variable, single term
x = Variable(randn(5))
b = randn(5)
cf = 10*norm(x-b,1)
xAll = RegLS.extract_variables(cf)
@test RegLS.is_proximable(cf) == true

f = RegLS.extract_proximable(xAll,cf)
@test norm(f(~x) - 10*norm(~x-b,1)) < 1e-12

# single variable, multiple terms with GetIndex
 x = Variable(randn(5))
 b = randn(2)
 cf = 10*norm(x[1:2]-b,1)+norm(x[3:5],2)
 xAll = RegLS.extract_variables(cf)
 @test RegLS.is_proximable(cf) == true
 f = RegLS.extract_proximable(xAll,cf)
 @test norm(f(~x) - sum([10*norm((~x)[1:2]-b,1);norm((~x)[3:5],2)])) < 1e-12

# multiple variables, multiple terms
x1 = Variable(randn(5))
b1 = randn(5)
x2 = Variable(randn(3))
b2 = randn(3)

cf = 10*norm(x2-b2,1)+norm(x1+b1,2)
xAll = (x1,x2)
@test RegLS.is_proximable(cf) == true
f = RegLS.extract_proximable(xAll,cf)
@test norm(f.fs[1](~x1)-norm(~x1+b1,2) ) < 1e-12
@test norm(f.fs[2](~x2)-10*norm(~x2-b2,1) ) < 1e-12

x1 = Variable(randn(5))
b1 = randn(5)
x2 = Variable(randn(5))
b2 = randn(5)

cf = 10*norm(x2+x1+b2,1)
xAll = (x1,x2)
@test RegLS.is_proximable(cf) == true
f = RegLS.extract_proximable(xAll,cf) #TODO

# multiple variables, missing terms
x1 = Variable(randn(5))
b1 = randn(5)
x2 = Variable(randn(3))
b2 = randn(3)

cf = 10*norm(x2-b2,1)
xAll = (x1,x2)
@test RegLS.is_proximable(cf) == true
f = RegLS.extract_proximable(xAll,cf)
@test f.fs[1](~x1) == 0.
@test norm(f.fs[2](~x2)-10*norm(~x2-b2,1) ) < 1e-12

# multiple variables, multiple terms, with GetIndex
x1 = Variable(randn(5))
b1 = randn(5)
x2 = Variable(randn(3))
b2 = randn(3)

cf = norm(x1[3:5]+b1[3:5],1)+10*norm(x2-b2,1)+norm(x1[1:2]+b1[1:2],2)
xAll = (x1,x2)
@test RegLS.is_proximable(cf) == true
f = RegLS.extract_proximable(xAll,cf)
@test norm(f.fs[1](~x1)-norm((~x1)[1:2]+b1[1:2],2)-norm((~x1)[3:5]+b1[3:5],1) ) < 1e-12
@test norm(f.fs[2](~x2)-10*norm(~x2-b2,1) ) < 1e-12

@printf("\nTesting problem construction primal\n")

x = Variable(10)
A = randn(5, 10)
y = Variable(7)
B = randn(5, 7)
b = randn(5)

prob = problem(ls(A*x + b), norm(x, 2) <= 1.0)
slv = RegLS.solve!(prob, PG())

~x .= 0.
prob = problem(ls(A*x - B*y + b) + norm(y, 1), norm(x, 2) <= 1.0)
slv = RegLS.solve!(prob, FPG())

@printf("\n Testing @minimize \n")
~x .= 0.
~y .= 0.
slv = @minimize ls(A*x - B*y + b) st norm(x, 2) <= 1e4, norm(y, 1) <= 1.0 with PG()
~x .= 0.
slv = @minimize ls(A*x - b) st norm(x, 1) <= 1.0 with PG()
~x .= 0.
slv = @minimize ls(A*x - b) st norm(x, 1) <= 1.0
~x .= 0.
slv = @minimize ls(A*x - b) + norm(x, 1) with PG()
~x .= 0.
slv = @minimize ls(A*x - b) + norm(x, 1)
~x .= 0.
slv = @minimize ls(A*x - b) 

x = Variable(5)
A = randn(10, 5)
b = randn(10)

#TODO many many more tests
@printf("\n Testing @minimize nonlinear \n")
slv = @minimize ls(sigmoid(A*x,10) - b)+norm(x,1) with PG() 
xpg = copy(~x)
~x .= 0.
slv = @minimize ls(sigmoid(A*x,10) - b)+norm(x,1) with ZeroFPR() 
xz = copy(~x)
~x .= 0.

@test norm(xz-xpg) <1e-7




