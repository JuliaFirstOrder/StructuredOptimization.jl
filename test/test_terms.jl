@printf("\nTesting cost terms\n")

# Simple Terms

x = Variable(10)
cf = norm(x, 0)
cf = norm(x, 0) <= 3
cf = norm(x, 1)
cf = norm(x, 1) <= 1.5
cf = norm(x, 2)
cf = norm(x, 2) <= 2.3
cf = norm(x, Inf)
cf = x <= 3.0
cf = 3.0 <= x
cf = x >= 1.0
cf = 1.0 >= x
cf = x in [-5.0, 5.0]
cf = norm(x)^2
cf = norm(x, 2)^2

X = Variable(10, 10)
cf = rank(X) <= 3

# Summing terms

x = Variable(10)
cf = ls(x) + norm(x, 1)

# More complex situations

x = Variable(10)
A = randn(5, 10)
y = Variable(7)
B = randn(5, 7)
b = randn(5)

cf1 = ls(A*x - b) + norm(x, 1)
cf2 = ls(A*x - B*y + b) + norm(y, 1) + norm(y, 2)
cf3 = 0.5*norm(A*x - B*y + b, 2)^2 + norm(x, 1) + norm(y, 2)

# u = Variable(5)
# w = Variable(5)
# z = Variable(5)
#
# cf4 = norm(A*x + z)
# @test RegLS.is_proximable(cf4) == false
# @test RegLS.is_smooth(cf4) == false
# @test RegLS.is_smooth(cf4^2) == true
#
# cf5 = norm(w + z)^2
# @test RegLS.is_proximable(cf5) == true
# @test RegLS.is_smooth(cf5) == true
#
# cf6 = norm(x, 1) + norm(y, 2)
# @test RegLS.is_proximable(cf6...) == true
# @test RegLS.is_smooth(cf6...) == false
#
# cf7 = norm(x, 1) + (x in [-1.0, +1.0])
# @test RegLS.is_proximable(cf7...) == false
# @test RegLS.is_smooth(cf7...) == false
#
# cf8 = norm(w + z, 1) + (norm(u) <= 1.0)
# @test RegLS.is_proximable(cf8...) == true
# @test RegLS.is_smooth(cf8...) == false

# testing tidy_up



A = randn(7,3)
n1,n2,n3,n4,n5 = 3,3,4,4,7
x1,x2,x3,x4,x5 = Variable(randn(n1)),Variable(randn(n2)),Variable(randn(n3)),Variable(randn(n4)),Variable(randn(n5))

cf = ls(x1+x2)
xAll = RegLS.get_all_variables(cf)
@test xAll == (x1,x2)

cf = ls(x1+x2)+ls(x1)
xAll = RegLS.get_all_variables(cf)
@test xAll == (x1,x2)

cf = ls(x1+x2)+ls(x3+x4)+ls(x5)+ls(x5+A*x2)+ls(x1)+ls(x5)
xAll = RegLS.get_all_variables(cf)
@test xAll == (x1,x2,x3,x4,x5)












