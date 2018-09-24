println("\nTesting linear expressions\n")

#### * ####
n, m1, m2, k = 3, 4, 5, 6
x1 = Variable(randn(m1))
x2 = Variable(randn(m2))
A1 = randn(n, m1)
A2 = randn(n, m2)
b  = randn(n)
b0  = pi
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
# multiply with Variable
ex1 = opA1*x1 
@test variables(ex1) == (x1,)
@test operator(ex1)*(~x1) == A1*(~x1)
B1 = randn(k, n)
opB1 = MatrixOp(B1)
# multiply with Expression
ex2 = opB1*ex1
@test variables(ex2) == (x1,)
@test norm(operator(ex2)*(~x1) - B1*A1*(~x1)) < 1e-12
# multiply with Expression with multiple variables
ex3 = opB1*(opA1*x1+opA2*x2)
@test variables(ex3) == (x1,x2)
@test norm(operator(ex3)*(~x1,~x2) - B1*(A1*(~x1)+A2*(~x2))) < 1e-12
# multiply with displacemented Array Expression with multiple variables
ex3 = opB1*(opA1*x1+opA2*x2+b)
@test variables(ex3) == (x1,x2)
@test norm(displacement(ex3) - B1*b) < 1e-12
# multiply with displacemented scalar Expression with multiple variables
ex3 = opB1*(opA1*x1+opA2*x2+b0)
@test variables(ex3) == (x1,x2)
@test norm(displacement(ex3) - B1*(ones(size(B1,2))*b0)) < 1e-12
@test_throws ArgumentError MatrixOp(randn(n,m1+1))*x1
@test_throws ArgumentError MatrixOp(randn(n,m1))*Variable(Complex{Float64},m1)

n, m1, m2, k = 3, 4, 5, 6
A1 = randn(n, m1)
A2 = randn(n, m2)
b1 = randn(n,n)
b2 = randn(n,n)
opA1 = MatrixOp(A1,n)
opA2 = MatrixOp(A2,n)
x1, x2  = Variable(randn(m1,n)), Variable(randn(m2,n))
# multiply Expressions (NonLinearCompose) 
ex = (opA1*x1)*(opA2*x2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1))*(A2*(~x2))) < 1e-12
ex = (opA1*x1-b1)*(opA2*x2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1)-b1)*(A2*(~x2))) < 1e-12
ex = (opA1*x1)*(opA2*x2+b2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1))*(A2*(~x2)+b2)) < 1e-12
ex = (opA1*x1+b1)*(opA2*x2+b2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1)+b1)*(A2*(~x2)+b2)) < 1e-12
@test_throws ErrorException (opA1*x1)*(opA1*x1)

n, m1, m2, k = 3, 4, 5, 6
A1 = randn(n, m1)
A2 = randn(n, m2)
b1 = randn(n,n)
b2 = randn(n,n)
opA1 = MatrixOp(A1,n)
opA2 = MatrixOp(A2,n)
x1, x2  = Variable(randn(m1,n)), Variable(randn(m2,n))
## multiply Expressions elementwise (Hadamard) 
ex = (opA1*x1).*(opA2*x2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1)).*(A2*(~x2))) < 1e-12
ex = (opA1*x1-b1).*(opA2*x2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1)-b1).*(A2*(~x2))) < 1e-12
ex = (opA1*x1).*(opA2*x2+b2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1)).*(A2*(~x2)+b2)) < 1e-12
ex = (opA1*x1+b1).*(opA2*x2+b2)
@test variables(ex) == (x1,x2)
@test norm(affine(ex)*(~x1,~x2) - (A1*(~x1)+b1).*(A2*(~x2)+b2)) < 1e-12
@test_throws ErrorException (opA1*x1)*(opA1*x1)

##### reshape ####
m,n = 8,10
A = randn(n,m)
x = Variable(randn(m))
b = randn(n)
B = reshape(b,2,5)

ex = reshape(x,4,2)
@test norm(operator(ex)*~x - reshape(~x,4,2)) < 1e-12

ex2 = reshape(A*x,2,5)
@test norm(operator(ex2)*~x - reshape(A*~x,2,5)) < 1e-12

ex3 = reshape(A*x,2,5)+B
@test norm(operator(ex2)*~x+displacement(ex3)- reshape(A*~x,2,5)-B) < 1e-12

ex4 = reshape(A*x-b,2,5)
@test norm(operator(ex4)*~x+displacement(ex4)- reshape(A*~x-b,2,5)) < 1e-12

#### + ####

# sum same variable
n, m = 3,4
x = Variable(randn(m))
A = randn(n, m)
B = randn(n, m)
opA = MatrixOp(A)
opB = MatrixOp(B)
ex1 = opA*x+opB*x 
@test variables(ex1) == (x,)
@test norm(operator(ex1)*~x - (opA+opB)*~x) < 1e-12

# sum different variables
n, m1, m2 = 3, 4, 5
xa = Variable(randn(m1))
xb = Variable(randn(m2))
xc = Variable(randn(n))
xd = Variable(randn(n))
A = randn(n, m1)
A2 = randn(n, m1)
B = randn(n, m2)
b = randn(n)
b0 = pi
opA = MatrixOp(A)
opA2 = MatrixOp(A2)
opB = MatrixOp(B)
opI = Eye(n)

# (+) sum different variables no HCAT
ex1 = opA*xa+opB*xb 
@test variables(ex1) == (xa,xb)
@test norm(operator(ex1)*(~xa,~xb) - hcat(opA,opB)*(~xa,~xb)) <1e-12

# (+) sum of same variables
ex2 = opA2*xa 
exs1 = ex1+ex2
exs2 = ex2+ex1
@test variables(exs1) == (xa,xb)
@test norm(operator(exs1)*(~xa,~xb) - hcat(opA+opA2,opB)*(~xa,~xb)) <1e-12
@test variables(exs2) == (xa,xb)
@test norm(operator(exs2)*(~xa,~xb) - hcat(opA+opA2,opB)*(~xa,~xb)) <1e-12

# (+) sum of different variables with HCAT
exs3 = exs1+exs2
@test variables(exs3) == (xa,xb)
@test norm(operator(exs3)*(~xa,~xb) - hcat(2*(opA+opA2),2*opB)*(~xa,~xb)) <1e-12
# (+) sum of different variables with HCAT
exs4 = exs1+(xc+xd)
@test variables(exs4) == (xa,xb,xc,xd)
@test norm(operator(exs4)*(~xa,~xb,~xc,~xd) - hcat(opA+opA2,opB,opI,opI)*(~xa,~xb,~xc,~xd)) <1e-12

# (+) sum Array
ex1 = xd+b
@test norm(displacement(ex1) - b) == 0.

# (+) sum scalar
ex2 = opB*xb+b0
@test (displacement(ex2) - b0) == 0.

# (+) sum displacemented expressions
ex3 = ex1+ex2
@test norm(displacement(ex3) - (b.+b0)) == 0.

#### (.+) sum 
n = 3
b = randn(n)

x1 = Variable(randn(1))
x2 = Variable(randn(n))
ex1 = x1.+x2 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9

x1 = Variable(randn(1))
x2 = Variable(randn(n))
ex1 = x1.+(x2+2) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9
@test displacement(ex1) == 2

x1 = Variable(randn(1))
x2 = Variable(randn(n))
ex1 = (x1+2).+(x2+b) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9
@test displacement(ex1) == (b.+2)

x1 = Variable(randn(n))
x2 = Variable(randn(1))
ex1 = x1.+x2 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9

x1 = Variable(randn(n))
x2 = Variable(randn(1))
ex1 = x1.+(x2+2) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9
@test displacement(ex1) == 2

x1 = Variable(randn(n))
x2 = Variable(randn(1))
ex1 = (x1+b).+(x2+2) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9
@test displacement(ex1) == (b.+2)

n,m =2,4
x1 = Variable(randn(n,m))
x2 = Variable(randn(1,m))
ex1 = x1.+x2+6
@test norm(operator(ex1)*(~x1,~x2)-((~x1).+(~x2))) < 1e-9
@test displacement(ex1) == 6

# #### (.-) sum 

n = 3
b = randn(n)

x1 = Variable(randn(1))
x2 = Variable(randn(n))
ex1 = x1.-x2 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9

x1 = Variable(randn(1))
x2 = Variable(randn(n))
ex1 = x1.-(x2+2) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9
@test displacement(ex1) == -2

x1 = Variable(randn(1))
x2 = Variable(randn(n))
ex1 = (x1+2).-(x2+b) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9
@test displacement(ex1) == (2 .-b)

x1 = Variable(randn(n))
x2 = Variable(randn(1))
ex1 = x1.-x2 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9

x1 = Variable(randn(n))
x2 = Variable(randn(1))
ex1 = x1.-(x2+2) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9
@test displacement(ex1) == -2

x1 = Variable(randn(n))
x2 = Variable(randn(1))
ex1 = (x1+b).-(x2+2) 
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9
@test displacement(ex1) == (b.-2)

n,m =2,4
x1 = Variable(randn(n,m))
x2 = Variable(randn(1,m))
ex1 = x1.-x2+6
@test norm(operator(ex1)*(~x1,~x2)-((~x1).-(~x2))) < 1e-9
@test displacement(ex1) == 6

# (-) sum different variables no HCAT
ex1 = opA*xa-opB*xb 
@test variables(ex1) == (xa,xb)
@test norm(operator(ex1)*(~xa,~xb) - hcat(opA,-opB)*(~xa,~xb)) <1e-12

# (-) sum of same variables
ex2 = opA2*xa 
exs1 = ex1-ex2
exs2 = ex2-ex1
@test variables(exs1) == (xa,xb)
@test norm(operator(exs1)*(~xa,~xb) - hcat(opA-opA2,-opB)*(~xa,~xb)) <1e-12
@test variables(exs2) == (xa,xb)
@test norm(operator(exs2)*(~xa,~xb) - hcat(-opA+opA2,+opB)*(~xa,~xb)) <1e-12

# (-) sum of same variables with HCAT
exs3 = exs1-exs2
@test variables(exs3) == (xa,xb)
@test norm(operator(exs3)*(~xa,~xb) - hcat(2*(opA-opA2),-2*opB)*(~xa,~xb)) < 1e-12
# (-) sum of different variables with HCAT
exs4 = exs1-(xc-xd)
@test variables(exs4) == (xa,xb,xc,xd)
@test norm(operator(exs4)*(~xa,~xb,~xc,~xd) - hcat(opA-opA2,-opB,-opI,opI)*(~xa,~xb,~xc,~xd)) <1e-12

# (-) sum Array
ex1 = xd-b
@test norm(displacement(ex1) + b) == 0.

# (-) sum scalar
ex2 = opB*xb-b0
@test (displacement(ex2) + b0) == 0.

# (+) sum displacemented expressions
ex3 = ex1-ex2
@test norm(displacement(ex3) - (-b.+b0)) == 0.

@test_throws DimensionMismatch MatrixOp(randn(10,20))*Variable(20)+randn(11)
@test_throws ErrorException MatrixOp(randn(10,20))*Variable(20)+(3+im)
