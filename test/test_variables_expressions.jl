@printf("\nTesting variables\n")

n, m, k = 3, 4, 5
x1 = Variable(n)
x1t = Variable(Float32, n)
x1i = Variable(randn(n))
x2 = Variable(n, m)
x2t = Variable(Float32, n, m)
x2i = Variable(randn(n, m))
x3 = Variable(n, m, k)
x3t = Variable(Float32, n, m, k)
x3i = Variable(randn(n, m, k))

@printf("\nTesting linear expressions\n")

A = randn(10, n)
opA = MatrixOp(A)
lex1 = LinearExpression(opA, x1)
x = randn(n)
@test vecnorm(lex1.L*x - A*x, Inf) <= 1e-12

B = randn(5, 10)
opB = MatrixOp(B)
lex2 = LinearExpression(opB, lex1)
@test vecnorm(lex2.L*x - B*(A*x), Inf) <= 1e-12

C = rand(Complex{Float64}, 10, n)
opC = MatrixOp(C)
x1c = Variable(Complex{Float64}, n)
lex3 = LinearExpression(opC, x1c)
x = rand(Complex{Float64}, n)
@test vecnorm(lex3.L*x - C*x, Inf) <= 1e-12

@printf("\nTesting affine expressions\n")

m, n1, n2 = 5, 7, 10
x = Variable(n1)
A = randn(m, n1)
lex1 = LinearExpression(MatrixOp(A), x)
aex1 = AffineExpression(lex1)
b = randn(m)
aex1 = AffineExpression(lex1, b)
y = Variable(n2)
B = randn(m, n2)
lex2 = LinearExpression(MatrixOp(B), y)

aex1 = lex1 + b
aex1 = lex1 - b
aex1 = b - lex1
aex1 = b + lex1
aex2 = aex1 + b
aex2 = aex1 - b
aex2 = b + aex1
aex2 = b - aex1
aex3 = aex1 + aex2
aex3 = aex1 - aex2
aex3 = aex2 + aex1
aex3 = aex2 - aex1
