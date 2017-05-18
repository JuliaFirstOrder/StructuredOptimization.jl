@printf("\nTesting linear operators calculus rules\n")

# test HCAT

m, n1, n2 = 4, 7, 5
A1 = randn(m, n1)
A2 = randn(m, n2)
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
opH = HCAT(opA1, opA2)
x1 = randn(n1)
x2 = randn(n2)
y1 = zeros(m)
A_mul_B!(y1, opH, (x1, x2))
y2 = A1*x1 + A2*x2

@test norm(y1-y2) <= 1e-12

# test VCAT of HCAT's

m = 5
A1 = randn(m, n1)
A2 = randn(m, n2)
opV = VCAT(opH, HCAT(MatrixOp(A1), MatrixOp(A2)))
y1 = (zeros(y1), zeros(m))
A_mul_B!(y1, opV, (x1, x2))
y2 = (y2, A1*x1 + A2*x2)

@test all(norm.(y1 .- y2) .<= 1e-12)

# test HCAT of VCAT's

n1, n2, m1, m2 = 3, 5, 4, 7
A = randn(m1, n1); opA = MatrixOp(A)
B = randn(m1, n2); opB = MatrixOp(B)
C = randn(m2, n1); opC = MatrixOp(C)
D = randn(m2, n2); opD = MatrixOp(D)
opV = HCAT(VCAT(opA, opC), VCAT(opB, opD))
x1 = randn(n1)
x2 = randn(n2)
y1 = (zeros(m1), zeros(m2))
A_mul_B!(y1, opV, (x1, x2))
y2 = (A*x1 + C*x2, B*x1 + D*x2)

@test all(norm.(y1 .- y2) .<= 1e-12)

# test Sum of HCAT's

m, n1, n2 = 4, 7, 5
A1 = randn(m, n1)
A2 = randn(m, n2)
B1 = randn(m, n1)
B2 = randn(m, n2)
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
opB1 = MatrixOp(B1)
opB2 = MatrixOp(B2)
opHA = HCAT(opA1, opA2)
opHB = HCAT(opB1, opB2)
opS = Sum(opHA, opHB)
x1 = randn(n1)
x2 = randn(n2)
y1 = zeros(m)
A_mul_B!(y1, opS, (x1, x2))
y2 = A1*x1 + B1*x1 + A2*x2 + B2*x2

@test norm(y1-y2) <= 1e-12

# test Sum of VCAT's

m1, m2, n = 4, 7, 5
A1 = randn(m1, n)
A2 = randn(m2, n)
B1 = randn(m1, n)
B2 = randn(m2, n)
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
opB1 = MatrixOp(B1)
opB2 = MatrixOp(B2)
opVA = VCAT(opA1, opA2)
opVB = VCAT(opB1, opB2)
opS = Sum(opVA, opVB)
x = randn(n)
y1 = (zeros(m1), zeros(m2))
A_mul_B!(y1, opS, x)
y2 = (A1*x + B1*x, A2*x + B2*x)

@test all(norm.(y1 .- y2) .<= 1e-12)
