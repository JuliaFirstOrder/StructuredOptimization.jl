using RegLS
using ProximalOperators
using Base.Test

srand(0)

m, n = 200, 100
A = randn(m, n)
b = randn(m)
var_x = Variable(n)
f = ls(A*var_x - b)
g = NormL1();
gamma = 1.0

q, s = RegLS.split_Quadratic(f)

x = deepcopy(~variable(f))

if ~isempty(q)
  res_q_x = RegLS.residual(q, x)
  res_q_xbar = deepcopy(res_q_x)
  grad_q_res, q_x = RegLS.gradient(q, res_q_x)
  grad_q_x = RegLS.At_mul_B(q, grad_q_res)
else
  grad_q_x = 0.0
  q_x = 0.0
end

if ~isempty(s)
  res_s_x = RegLS.residual(s, x)
  res_s_xbar = deepcopy(res_s_x)
  grad_s_res, s_x = RegLS.gradient(s, res_s_x)
  grad_s_x = RegLS.At_mul_B(q, grad_s_res)
else
  grad_s_x = 0.0
  s_x = 0.0
end

f_x = q_x + s_x
grad_f_x = grad_q_x + grad_s_x

gradstep_x = x - gamma*grad_f_x
xbar, g_xbar = prox(g, gradstep_x, gamma)

if ~isempty(q)
  RegLS.residual!(res_q_xbar, q, xbar)
  q_xbar = RegLS.cost(q, res_q_xbar)
else q_xbar = 0.0 end

if ~isempty(s)
  RegLS.residual!(res_s_xbar, s, xbar)
  s_xbar = RegLS.cost(s, res_s_xbar)
else s_xbar = 0.0 end

f_xbar = q_xbar + s_xbar
