println("\nTesting variables\n")

n, m, k = 3, 4, 5
x1t = Variable(Float32, n)
x1 = Variable(n)
x1i = Variable(randn(n))
x2 = Variable(n, m)
x2t = Variable(Float32, n, m)
x2i = Variable(randn(n, m))
x3 = Variable(n, m, k)
x3t = Variable(Float32, n, m, k)
xx = randn(n,m,k)
x3i = Variable(xx)

@test eltype(x3i) == eltype(xx)
@test size(x3i) == size(xx)
@test size(x3i,1) == size(xx,1)
@test xx == (~x3i)

@test typeof(operator(x1)) <: Eye
@test variables(x1) == x1
