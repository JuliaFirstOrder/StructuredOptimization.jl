using RegLS
using Base.Test

@printf("Testing L-BFGS routines\n")

Q = [32.0000 13.1000 -4.9000 -3.0000  6.0000  2.2000  2.6000  3.4000 -1.9000 -7.5000;
 13.1000 18.3000 -5.3000 -9.5000  3.0000  2.1000  3.9000  3.0000 -3.6000 -4.4000;
 -4.9000 -5.3000  7.7000  2.1000 -0.4000 -3.4000 -0.8000 -3.0000  5.3000  5.5000;
 -3.0000 -9.5000  2.1000 20.1000  1.1000  0.8000  -12.4000 -2.5000  5.5000  2.1000;
  6.0000  3.0000 -0.4000  1.1000  3.8000  0.6000  0.5000  0.9000 -0.4000 -2.0000;
  2.2000  2.1000 -3.4000  0.8000  0.6000  7.8000  2.9000 -1.3000 -4.3000 -5.1000;
  2.6000  3.9000 -0.8000  -12.4000  0.5000  2.9000 14.5000  1.7000 -4.9000  1.2000;
  3.4000  3.0000 -3.0000 -2.5000  0.9000 -1.3000  1.7000  6.6000 -0.8000  2.7000;
 -1.9000 -3.6000  5.3000  5.5000 -0.4000 -4.3000 -4.9000 -0.8000  7.9000  5.7000;
 -7.5000 -4.4000  5.5000  2.1000 -2.0000 -5.1000  1.2000  2.7000  5.7000 16.1000]

q = [2.9000, 0.8000, 1.3000, -1.1000, -0.5000, -0.3000, 1.0000, -0.3000, 0.7000, -2.1000]

xs = [1.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09;
      0.09 1.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08;
      0.08 0.09 1.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07;
      0.07 0.08 0.09 1.0 0.01 0.02 0.03 0.04 0.05 0.06;
      0.06 0.07 0.08 0.09 1.0 0.01 0.02 0.03 0.04 0.05]'

dirs_ref = [
  -3.476000000000000e+01    -6.861170733797231e-01    -1.621334774299757e-01    -2.008976150849174e-01    -2.317011191832649e-01;
  -1.367700000000000e+01    -1.661270665201917e+00     2.870743130038511e-01     2.237224648542354e-01     2.980080835636926e-02;
   2.961000000000000e+00     2.217225828759783e-01    -5.485761164147891e-01     4.811889625788801e-02    -1.267017945785352e-01;
   3.756000000000000e+00     5.615134140894827e-01     9.992734938824949e-02    -6.855884193567087e-01     4.328230970765587e-02;
  -5.618000000000001e+00    -1.922426760799171e-01    -1.332550298134261e-02    -2.729265954345345e-02    -2.437461022925742e-01;
  -1.571000000000000e+00    -8.961101045874649e-02     5.326252573648003e-02     3.651730112313705e-02     1.349716200511426e-02;
  -4.121000000000000e+00    -3.044802963260585e-01    -6.299408068289100e-02     6.325330777317102e-02    -7.155992987801297e-04;
  -3.709000000000000e+00    -1.996235459345302e-01     1.525398352758626e-02     2.871281112230844e-02    -3.513449694839536e-03;
   4.010000000000000e-01     1.267604425710271e-01    -7.776943954825602e-02    -1.285590864125103e-01    -5.603489763638488e-02;
   7.639999999999999e+00     3.360845247013288e-01    -2.335884953507600e-02    -3.204963735369062e-03     5.612114259243499e-02
]

dirs = zeros(10, 5) # matrix of directions (to be filled in)

mem = 3;
col = 0; # last column of Sk, Yk that was filled in
currmem = 0;
H0 = 1.0

lbfgs = RegLS.LBFGS.create(mem,Array{Float64,1})
x_old = 0;
grad_old = 0;

for i=1:5
    x = xs[:,i]
    grad = Q*x + q
    if i > 1
        Sk = x-x_old
        Yk = grad-grad_old
        RegLS.LBFGS.push(lbfgs, Sk, Yk)
        H0 = dot(Sk[:],Yk[:])/dot(Yk[:],Yk[:])
    end
    dir = RegLS.LBFGS.matvec(lbfgs, H0, -grad)
    dirs[:, i] = dir
    x_old = x;
    grad_old = grad;
end

@test vecnorm(dirs-dirs_ref, Inf)/(1+vecnorm(dirs_ref, Inf)) <= 1e-15
