@printf("\nTesting linear operators calculus rules\n")

verb = true

# test HCAT

m, n1, n2 = 4, 7, 5
A1 = randn(m, n1)
A2 = randn(m, n2)
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
opH = HCAT(opA1, opA2)
x1 = randn(n1)
x2 = randn(n2)
y1 = test_op(opH, (x1, x2), randn(m), verb)
y2 = A1*x1 + A2*x2
@test vecnorm(y1-y2) <= 1e-12

# test VCAT of HCAT's

m = 5
A1 = randn(m, n1)
A2 = randn(m, n2)
opV = VCAT(opH, HCAT(MatrixOp(A1), MatrixOp(A2)))
y1 = test_op(opV, (x1, x2), (randn(size(y1)), randn(m)), verb)
y2 = (y2, A1*x1 + A2*x2)

@test all(vecnorm.(y1 .- y2) .<= 1e-12)

# test HCAT of VCAT's

n1, n2, m1, m2 = 3, 5, 4, 7
A = randn(m1, n1); opA = MatrixOp(A)
B = randn(m1, n2); opB = MatrixOp(B)
C = randn(m2, n1); opC = MatrixOp(C)
D = randn(m2, n2); opD = MatrixOp(D)
opV = HCAT(VCAT(opA, opC), VCAT(opB, opD))
x1 = randn(n1)
x2 = randn(n2)
y1 = test_op(opV, (x1, x2), (randn(m1), randn(m2)), verb)
y2 = (A*x1 + B*x2, C*x1 + D*x2)

@test all(vecnorm.(y1 .- y2) .<= 1e-12)

# test Sum of HCAT's

m, n1, n2, n3 = 4, 7, 5, 3
A1 = randn(m, n1)
A2 = randn(m, n2)
A3 = randn(m, n3)
B1 = randn(m, n1)
B2 = randn(m, n2)
B3 = randn(m, n3)
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
opA3 = MatrixOp(A3)
opB1 = MatrixOp(B1)
opB2 = MatrixOp(B2)
opB3 = MatrixOp(B3)
opHA = HCAT(opA1, opA2, opA3)
opHB = HCAT(opB1, opB2, opB3)
opS = Sum(opHA, opHB)
x1 = randn(n1)
x2 = randn(n2)
x3 = randn(n3)
y1 = test_op(opS, (x1, x2, x3), randn(m), verb)
y2 = A1*x1 + B1*x1 + A2*x2 + B2*x2 + A3*x3 + B3*x3

@test vecnorm(y1-y2) <= 1e-12

# test Sum of VCAT's

m1, m2, n = 4, 7, 5
A1 = randn(m1, n)
A2 = randn(m2, n)
B1 = randn(m1, n)
B2 = randn(m2, n)
C1 = randn(m1, n)
C2 = randn(m2, n)
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
opB1 = MatrixOp(B1)
opB2 = MatrixOp(B2)
opC1 = MatrixOp(C1)
opC2 = MatrixOp(C2)
opVA = VCAT(opA1, opA2)
opVB = VCAT(opB1, opB2)
opVC = VCAT(opC1, opC2)
opS = Sum(opVA, opVB, opVC)
x = randn(n)
y1 = test_op(opS, x, (randn(m1), randn(m2)), verb)
y2 = (A1*x + B1*x +C1*x, A2*x + B2*x + C2*x)

@test all(vecnorm.(y1 .- y2) .<= 1e-12)
