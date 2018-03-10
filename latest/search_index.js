var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#StructuredOptimization.jl-1",
    "page": "Home",
    "title": "StructuredOptimization.jl",
    "category": "section",
    "text": "StructuredOptimization.jl is a high-level modeling language that utilizes a syntax that is very close to the mathematical formulation of an optimization problem.This user-friendly interface acts as a parser to utilize three different packages:ProximalOperators.jl provides proximal mappings of functions that are frequently used in signal processing and optimization.\nAbstractOperators.jl provides algorithms for the evaluation and combination of forward and (Jacobian) adjoint of linear and nonlinear mappings.\nProximalAlgorithms.jl is a library of proximal algorithms (aka splitting algorithms) solvers.StructuredOptimization.jl can handle large-scale convex and nonconvex problems with nonsmooth cost functions. It supports complex variables as well. See the Quick tutorial guide and the Demos."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "From the Julia command line hit Pkg.clone(\"https://github.com/kul-forbes/StructuredOptimization.jl.git\"). Once the package is installed you can update it along with the others issuing Pkg.update() in the command line."
},

{
    "location": "index.html#Citing-1",
    "page": "Home",
    "title": "Citing",
    "category": "section",
    "text": "If you use StructuredOptimization.jl for published work, we encourage you to cite:N. Antonello, L. Stella, P. Patrinos, T. van Waterschoot, “Proximal Gradient Algorithms: Applications in Signal Processing,” arXiv:1803.01621 (2018)."
},

{
    "location": "index.html#Credits-1",
    "page": "Home",
    "title": "Credits",
    "category": "section",
    "text": "StructuredOptimization.jl is developed by Lorenzo Stella and Niccolò Antonello at KU Leuven, ESAT/Stadius."
},

{
    "location": "tutorial.html#",
    "page": "Quick Tutorial Guide",
    "title": "Quick Tutorial Guide",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial.html#Quick-tutorial-guide-1",
    "page": "Quick Tutorial Guide",
    "title": "Quick tutorial guide",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial.html#Standard-problem-formulation-1",
    "page": "Quick Tutorial Guide",
    "title": "Standard problem formulation",
    "category": "section",
    "text": "Currently with StructuredOptimization.jl you can solve problems of the formunderset mathbfx textminimize  f(mathbfx) + g(mathbfx)where f is a smooth function while g is possibly nonsmooth."
},

{
    "location": "tutorial.html#Unconstrained-optimization-1",
    "page": "Quick Tutorial Guide",
    "title": "Unconstrained optimization",
    "category": "section",
    "text": "The least absolute shrinkage and selection operator (LASSO) belongs to this class of problems:underset mathbfx textminimize  tfrac12  mathbfA mathbfx - mathbfy ^2+ lambda  mathbfx _1Here the squared norm tfrac12  mathbfA mathbfx - mathbfy ^2 is a smooth function f wherelse the l_1-norm is a nonsmooth function g.This problem can be solved with only few lines of code:julia> using StructuredOptimization\n\njulia> n, m = 100, 10;                # define problem size\n\njulia> A, y = randn(m,n), randn(m);   # random problem data\n\njulia> x = Variable(n);               # initialize optimization variable\n\njulia> λ = 1e-2*norm(A\'*y,Inf);       # define λ    \n\njulia> @minimize ls( A*x - y ) + λ*norm(x, 1); # solve problem\n\njulia> ~x                             # inspect solution\n100-element Array{Float64,1}:\n  0.0\n  0.0\n  0.0\n  0.440254\n  0.0\n  0.0\n  0.0\n[...]note: Note\nThe function ls is a short hand notation for 0.5*norm(...)^2, namely a least squares term.It is possible to access to the solution by typing ~x.By default variables are initialized by Arrays of zeros.Different initializations can be set during construction x = Variable( [1.; 0.; ...] ) or by assignement ~x .= [1.; 0.; ...]."
},

{
    "location": "tutorial.html#Constrained-optimization-1",
    "page": "Quick Tutorial Guide",
    "title": "Constrained optimization",
    "category": "section",
    "text": "Constrained optimization is also encompassed by the Standard problem formulation:for a nonempty set mathcalS the constraint ofbeginalign*\nunderset mathbfx textminimize    f(mathbfx) \ntextsubject to   mathbfx in mathcalS\nendalign*can be converted into an indicator functiong(mathbfx) = delta_mathcalS (mathbfx) =  begincases\n    0        textif  mathbfx in mathcalS\n    +infty  textotherwise\n    endcasesto obtain the standard form. Constraints are treated as nonsmooth functions.This conversion is automatically performed by StructuredOptimization.jl.For example, the non-negative deconvolution problem:beginalign*\nunderset mathbfx textminimize    tfrac12  mathbfx * mathbfh - mathbfy ^2 \ntextsubject to   mathbfx geq 0\nendalign*where * stands fof convoluton and mathbfh contains the taps of a finite impluse response, can be solved using the following lines of code:julia> n = 10;                        # define problem size\n\njulia> x = Variable(n);               # define variable\n\njulia> h, y = randn(n), randn(2*n-1); # random filter taps and output\n\njulia> @minimize ls(conv(x,h)-y) st x >= 0.\nnote: Note\nThe convolution mapping was applied to the variable x using conv.StructuredOptimization.jl provides a set of functions that can be used to apply specific operators to variables and create mathematical expression.The available functions can be found in Mappings.In general it is more convenient to use these functions instead of matrices, as these functions apply efficient algorithms for the forward and adjoint mappings leading to matrix free optimization."
},

{
    "location": "tutorial.html#Using-multiple-variables-1",
    "page": "Quick Tutorial Guide",
    "title": "Using multiple variables",
    "category": "section",
    "text": "It is possible to use multiple variables which are allowed to be matrices or even tensors.For example a non-negative matrix factorization problem:beginalign*\nunderset mathbfX_1 mathbfX_2  textminimize    tfrac12  mathbfX_1 mathbfX_2 - mathbfY  \ntextsubject to   mathbfX_1 geq 0   mathbfX_2 geq 0\nendalign*can be solved using the following code:# matrix variables initialized with random coefficients\njulia> X1, X2 = Variable(rand(n,l)), Variable(rand(l,m));\n\njulia> Y = rand(n,m);\n\njulia> @minimize ls(X1*X2-Y) st X1 >= 0., X2 >= 0.\n"
},

{
    "location": "tutorial.html#Limitations-1",
    "page": "Quick Tutorial Guide",
    "title": "Limitations",
    "category": "section",
    "text": "Currently StructuredOptimization.jl supports only Proximal Gradient (aka Forward Backward) algorithms, which require specific properties of the nonsmooth functions and costraint to be applicable.In particular, the nonsmooth functions must lead to an efficiently computable proximal mapping.If we express the nonsmooth function g as the composition of a function tildeg with a linear operator A:g(mathbfx) =\ntildeg(A mathbfx)then a proximal mapping of g is efficiently computable if it satisifies the following properties:the mapping A must be a tight frame  namely it must satisfy A A^* = mu Id, where mu geq 0 and A^* is the adjoint of A and Id is the identity operator.\nif A is not a tight frame, than it must be possible write g as a separable sum g(mathbfx) =  sum_j h_j (B_j mathbfx_j) with mathbfx_j being a non-overlapping slices of mathbfx and B_j being tight frames.Let us analyze these rules with a series of examples.The LASSO example above satisfy the first rule:julia> @minimize ls( A*x - y ) + λ*norm(x, 1)\nsince the non-smooth function lambda  cdot _1 is not composed with any operator (or equivalently is composed with Id which is a tight frame).Also the following problem would be accepted:julia> @minimize ls( A*x - y ) + λ*norm(dct(x), 1)\nsince the discrete cosine transform (DCT) is orthogonal and is therefore a tight frame.On the other hand, the following problemjulia> @minimize ls( A*x - y ) + λ*norm(x, 1) st x >= 1.0\ncannot be solved through proximal gradient algorithms, since the second rule would be violated. Here the constraint would be converted into an indicator function and the nonsmooth function g can be written as the sum:g(mathbfx) =lambda  mathbfx _1 + delta_mathcalS (mathbfx)which is not separable.On the other hand this problem would be accepted:julia> @minimize ls( A*x - y ) + λ*norm(x[1:div(n,2)], 1) st x[div(n,2)+1:n] >= 1.0\nas not the optimization variables mathbfx are partitioned into non-overlapping groups.note: Note\nWhen the problem is not accepted it might be still possible to solve it: see Smoothing and Duality."
},

{
    "location": "expressions.html#",
    "page": "Expressions",
    "title": "Expressions",
    "category": "page",
    "text": ""
},

{
    "location": "expressions.html#Expressions-1",
    "page": "Expressions",
    "title": "Expressions",
    "category": "section",
    "text": "With StructuredOptimization.jl you can easily create mathematical expressions. Firstly, Variables must be defined: various Mappings can then be applied  following the application of Functions and constraints to create the Terms  that define the optimization problem. "
},

{
    "location": "expressions.html#Variables-1",
    "page": "Expressions",
    "title": "Variables",
    "category": "section",
    "text": ""
},

{
    "location": "expressions.html#StructuredOptimization.Variable",
    "page": "Expressions",
    "title": "StructuredOptimization.Variable",
    "category": "type",
    "text": "Variable([T::Type,] dims...)\n\nReturns a Variable of dimension dims initialized with an array of all zeros.\n\nVariable(x::AbstractArray)\n\nReturns a Variable of dimension size(x) initialized with x\n\n\n\n"
},

{
    "location": "expressions.html#Creating-Variables-1",
    "page": "Expressions",
    "title": "Creating Variables",
    "category": "section",
    "text": "Variablenote: Note\nStructuredOptimization.jl supports complex variables. It is possible to create them by specifying the type  Variable(Complex{Float64}, 10) or by initializing them with a complex array Variable(randn(10)+im*randn(10))."
},

{
    "location": "expressions.html#Base.:~",
    "page": "Expressions",
    "title": "Base.:~",
    "category": "function",
    "text": "~(x::Variable)\n\nReturns the Array of the variable x\n\n\n\n"
},

{
    "location": "expressions.html#Base.size",
    "page": "Expressions",
    "title": "Base.size",
    "category": "function",
    "text": "size(x::Variable, [dim...])\n\nLike size(A::AbstractArray, [dims...]) returns the tuple containing the dimensions of the variable x.\n\n\n\n"
},

{
    "location": "expressions.html#Base.eltype",
    "page": "Expressions",
    "title": "Base.eltype",
    "category": "function",
    "text": "eltype(x::Variable)\n\nLike eltype(x::AbstractArray) returns the type of the elements of x.\n\n\n\n"
},

{
    "location": "expressions.html#Utilities-1",
    "page": "Expressions",
    "title": "Utilities",
    "category": "section",
    "text": "~\nsize\neltype"
},

{
    "location": "expressions.html#Base.:+",
    "page": "Expressions",
    "title": "Base.:+",
    "category": "function",
    "text": "+(ex1::AbstractExpression, ex2::AbstractExpression)\n\nAdd two expressions. \n\nExamples\n\njulia> x,y = Variable(5), Variable(5)\n(Variable(Float64, (5,)), Variable(Float64, (5,)))\n\njulia> ex1 = x+y\n\njulia> z = Variable(2)\nVariable(Float64, (2,))\n\njulia> ex2 = randn(5,2)*z\n\n\nNotice that in order for two expressions to be added toghether their associated AbstractOperator  must have the same codomain:\n\njulia> operator(ex1)\n[I,I]  ℝ^5  ℝ^5 -> ℝ^5 \n\njulia> operator(ex2)\n▒  ℝ^2 -> ℝ^5 \n\njulia> ex3 = ex1 + ex2\n\n\nIt is also possible to use broadcasted addition:\n\njulia> z = Variable(1)\nVariable(Float64, (1,))\n\njulia> ex3.+z\n\n\n\n\n+(ex::AbstractExpression, b::Union{AbstractArray,Number})\n\nAdd a scalar or an Array to an expression: \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = x+4\n\n\nNotice that in order to add an array to ex, b must belong to the codomain  of the associated AbstractOperator of ex. \n\njulia> b = randn(10);\n\njulia> size(b), eltype(b)\n((10,), Float64)\n\njulia> size(operator(ex),1), codomainType(operator(ex))\n((10,), Float64)\n\njulia> ex + b\n\n\n\n\n"
},

{
    "location": "expressions.html#Summing-expressions-1",
    "page": "Expressions",
    "title": "Summing expressions",
    "category": "section",
    "text": "+"
},

{
    "location": "expressions.html#Base.:*",
    "page": "Expressions",
    "title": "Base.:*",
    "category": "function",
    "text": "*(A::AbstractOperator, ex::AbstractExpression)\n\nMultiply an \'AbstractExpressionby anAbstractOperator`. \n\nExample\n\njulia> A = FiniteDiff(Float64, (10,))\nδx  ℝ^10 -> ℝ^9\n\njulia> x = Variable(10);\n\njulia> ex = A*x;\n\njulia> B = DCT(9)\nℱc  ℝ^9 -> ℝ^9\n\njulia> ex2 = B*ex;\n\njulia> operator(ex2)\nℱc*δx  ℝ^10 -> ℝ^9\n\n\n\n*(A::AbstractMatrix, ex::AbstractExpression)\n\nMultiply an AbstractExpression by an AbstractMatrix. \n\njulia> A = randn(10,5);\n\njulia> x = Variable(5)\nVariable(Float64, (5,))\n\njulia> A*x\n\n\nOther types of multiplications are also possible:\n\nleft array multiplication\njulia> X = Variable(10,5)\nVariable(Float64, (10, 5))\n\njulia> X*randn(5,6)\n\nscalar multiplication:\njulia> π*X\n\nelementwise multiplication:\njulia> randn(10,5).*X\n\n\n\n\n*(A::AbstractExpression, ex::AbstractExpression)\n\nMultiply an AbstractExpression by another AbstractExpression. \n\nExamples\n\njulia> W1 = Variable(10,5)\nVariable(Float64, (10, 5))\n\njulia> W2 = Variable(5,15)\nVariable(Float64, (5, 15))\n\njulia> ex = W1*σ(W2);\n\njulia> operator(ex)\nI*σ  ℝ^(10, 5)  ℝ^(5, 15) -> ℝ^(10, 15)\n\n\n\n\n"
},

{
    "location": "expressions.html#Multiplying-expressions-1",
    "page": "Expressions",
    "title": "Multiplying expressions",
    "category": "section",
    "text": "*"
},

{
    "location": "expressions.html#Mappings-1",
    "page": "Expressions",
    "title": "Mappings",
    "category": "section",
    "text": "As shown in the Quick tutorial guide it is possible to apply different mappings to the variables  using a simple syntax. Alternatively, as shown in Multiplying expressions, it is possible to define the mappings using  AbstractOperators.jl and to apply them  to the variable (or expression) through multiplication."
},

{
    "location": "expressions.html#Base.getindex",
    "page": "Expressions",
    "title": "Base.getindex",
    "category": "function",
    "text": "getindex(x::AbstractExpression, dims...)\n\nSlices the codomain of ex.\n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> getindex(x,1:4)\n\njulia> X = Variable(10,4)\nVariable(Float64, (10, 4))\n\njulia> X[:,1:3]\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.reshape",
    "page": "Expressions",
    "title": "Base.reshape",
    "category": "function",
    "text": "reshape(x::AbstractExpression, dims...)\n\nReshapes the codomain of the expression.\n\nExample\n\njulia> A,b = randn(10,3), randn(10);\n\njulia> reshape(A*x-b,2,5)\n\n\n\n\n"
},

{
    "location": "expressions.html#Basic-mappings-1",
    "page": "Expressions",
    "title": "Basic mappings",
    "category": "section",
    "text": "getindex\nreshape"
},

{
    "location": "expressions.html#Base.DFT.fft",
    "page": "Expressions",
    "title": "Base.DFT.fft",
    "category": "function",
    "text": "fft(x::AbstractExpression)\n\nApplies the Fourier transform. See documentation for Base.fft and AbstractOperator.DFT.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = fft(x)\n\njulia> operator(ex)\nℱ  ℝ^10 -> ℂ^10\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DFT.ifft",
    "page": "Expressions",
    "title": "Base.DFT.ifft",
    "category": "function",
    "text": "ifft(x::AbstractExpression)\n\nApplies the inverse Fourier transform. See documentation for Base.ifft and AbstractOperator.IDFT.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = ifft(x)\n\njulia> operator(ex)\nℱ⁻¹  ℝ^10 -> ℂ^10\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DFT.rfft",
    "page": "Expressions",
    "title": "Base.DFT.rfft",
    "category": "function",
    "text": "rfft(x::AbstractExpression, [, dims] )\n\nApplies the Fourier transform exploiting the conjugate symmetry.  See documentation for Base.rfft and AbstractOperator.RDFT.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = rfft(x)\n\njulia> operator(ex)\nℱ  ℝ^10 -> ℂ^6\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DFT.irfft",
    "page": "Expressions",
    "title": "Base.DFT.irfft",
    "category": "function",
    "text": "irfft(x::AbstractExpression, d, [, dims] )\n\nApplies the inverse Fourier transform exploiting the conjugate symmetry. \n\nd must indicate the length of the real codomain.\n\nSee documentation for Base.irfft and AbstractOperator.IRDFT.  \n\nExample\n\njulia> x = Variable(Complex{Float64}, 10)\nVariable(Float64, (10,))\n\njulia> ex = irfft(x,19)\n\njulia> operator(ex)\nℱ⁻¹  ℂ^10 -> ℝ^19\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DFT.FFTW.dct",
    "page": "Expressions",
    "title": "Base.DFT.FFTW.dct",
    "category": "function",
    "text": "dct(x::AbstractExpression)\n\nApplies the discrete cosine transform. See documentation for Base.dct and AbstractOperator.DCT.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = dct(x)\n\njulia> operator(ex)\nℱc  ℝ^10 -> ℝ^10\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DFT.FFTW.idct",
    "page": "Expressions",
    "title": "Base.DFT.FFTW.idct",
    "category": "function",
    "text": "idct(x::AbstractExpression)\n\nApplies the inverse discrete cosine transform. See documentation for Base.idct and AbstractOperator.IDCT.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = idct(x)\n\njulia> operator(ex)\nℱc⁻¹  ℝ^10 -> ℝ^10\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DSP.conv",
    "page": "Expressions",
    "title": "Base.DSP.conv",
    "category": "function",
    "text": "conv(x::AbstractExpression, h::AbstractVector)\n\nPerform discrete convolution with h. See documentation for Base.conv and AbstractOperator.Conv.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = conv(x, randn(5))\n\njulia> operator(ex)\n★  ℝ^10 -> ℝ^14\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DSP.xcorr",
    "page": "Expressions",
    "title": "Base.DSP.xcorr",
    "category": "function",
    "text": "xcorr(x::AbstractExpression, h::AbstractVector)\n\nPerforms the cross correlation of the codomain of ex with h. See documentation for Base.xcorr and AbstractOperator.Xcorr.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = xcorr(x, randn(5))\n\njulia> operator(ex)\n◎  ℝ^10 -> ℝ^19\n\n\n\n\n"
},

{
    "location": "expressions.html#Base.DSP.filt",
    "page": "Expressions",
    "title": "Base.DSP.filt",
    "category": "function",
    "text": "filt(x::AbstractExpression, b::AbstractVector, [a::AbstractVector])\n\nFilter with the finite impulse response (FIR) b. Alternatively infinite impulse responses filters (IIR) can be used as well by specifying the coefficients a.\n\nSee documentation for Base.filt and AbstractOperator.Filt.  \n\nExample\n\njulia> x = Variable(10)\nVariable(Float64, (10,))\n\njulia> ex = filt(x, [1.;0.;0.],[1.;0.;1.])\n\njulia> operator(ex)\nIIR  ℝ^10 -> ℝ^10\n\njulia> ex = filt(x, [1.;0.;0.])\n\njulia> operator(ex)\nFIR  ℝ^10 -> ℝ^10\n\n\n\n\n"
},

{
    "location": "expressions.html#StructuredOptimization.mimofilt",
    "page": "Expressions",
    "title": "StructuredOptimization.mimofilt",
    "category": "function",
    "text": "mimofilt(X::AbstractExpression, B::Vector{AbstractVector}, [A::Vector{AbstractVector}])\n\nMultiple-input multiple-output filter: like filt but allows for multipule inputs. \n\nmathbfy_i = sum_j = 1^M mathbfh_ij * mathbfx_j \n\nwhere mathbfy_i and mathbfx_j are the i-th and j-th columns of the output Y and input X matrices respectively.\n\nThe filters mathbfh_ij can be represented either by providing coefficients B and A (IIR) or B alone (FIR). These coefficients must be given in a Vector of Vectors. \n\nFor example for a 3 by 2 MIMO system (i.e. size(X,2) == 3 inputs and size(Y,2) == 2 outputs) B must be:\n\nB = [b11, b12, b13, b21, b22, b23]\n\nwhere bij are vector containing the filter coeffients of h_{i,j}.\n\nSee documentation of AbstractOperator.MIMOFilt.  \n\n\n\n"
},

{
    "location": "expressions.html#StructuredOptimization.zeropad",
    "page": "Expressions",
    "title": "StructuredOptimization.zeropad",
    "category": "function",
    "text": "zeropad(x::AbstractExpression, zp::Tuple)\n\nPerforms zeropadding:\n\nzp(mathbfx) = mathbfx 0 dots 0  ^T\n\nwhere zp is a tuple containing the number of nonzero elements that are added in each dimension.\n\nSee documentation of AbstractOperator.ZeroPad.  \n\n\n\n"
},

{
    "location": "expressions.html#DSP-mappings-1",
    "page": "Expressions",
    "title": "DSP mappings",
    "category": "section",
    "text": "fft\nifft\nrfft\nirfft\ndct\nidct\nconv\nxcorr\nfilt\nmimofilt\nzeropad"
},

{
    "location": "expressions.html#StructuredOptimization.finitediff",
    "page": "Expressions",
    "title": "StructuredOptimization.finitediff",
    "category": "function",
    "text": "finitediff(X::AbstractExpression, dir = 1)\n\nPerforms the discretized gradient over the specified direction dir obtained using forward finite differences. \n\nExample\n\njulia> X = Variable(10,10)\nVariable(Float64, (10,10))\n\njulia> ex = finitediff(X)\n\njulia> operator(ex)\nδx  ℝ^(10, 10) -> ℝ^(9, 10)\n\njulia> ex = finitediff(X,2)\n\njulia> operator(ex)\nδy  ℝ^(10, 10) -> ℝ^(10, 9)\n\n\n\n\n"
},

{
    "location": "expressions.html#StructuredOptimization.variation",
    "page": "Expressions",
    "title": "StructuredOptimization.variation",
    "category": "function",
    "text": "variation(X::AbstractExpression)\n\nReturns the variation of X using forward finite differences.  Specifically if X is in ℝ^(n, m) the codomain of variation(X) will  consist of ℝ^(n*m,2), having in the i-th column the vectorized gradient over the i-th direction. \n\nSee documentation of AbstractOperator.Variation.  \n\nExample\n\njulia> X = Variable(4,5,6)\nVariable(Float64, (4,5,6))\n\njulia> ex = variation(X)\n\njulia> operator(ex)\nƲ  ℝ^(4, 5, 6) -> ℝ^(120, 3)\n\n\n\n\n"
},

{
    "location": "expressions.html#Finite-differences-mappings-1",
    "page": "Expressions",
    "title": "Finite differences mappings",
    "category": "section",
    "text": "finitediff\nvariation"
},

{
    "location": "expressions.html#StructuredOptimization.sigmoid",
    "page": "Expressions",
    "title": "StructuredOptimization.sigmoid",
    "category": "function",
    "text": "sigmoid(x::AbstractExpression, γ = 1.0) σ(x::AbstractExpression, γ = 1.0)\n\nSigmoid function:\n\nsigma(mathbfx) = frac11+e^-gamma mathbfx  \n\nSee documentation of AbstractOperator.Sigmoid.  \n\n\n\n"
},

{
    "location": "expressions.html#Nonlinear-mappings-1",
    "page": "Expressions",
    "title": "Nonlinear mappings",
    "category": "section",
    "text": "sigmoid"
},

{
    "location": "expressions.html#StructuredOptimization.variables",
    "page": "Expressions",
    "title": "StructuredOptimization.variables",
    "category": "function",
    "text": "variables(ex::Expression)\n\nReturns a tuple containing the Variables of expression ex.\n\nExample\n\njulia> x,y = Variable(2),Variable(2);\n\njulia> ex = x+y;\n\njulia> variables(ex)\n(Variable(Float64, (2,)), Variable(Float64, (2,)))\n\n\n\n\n"
},

{
    "location": "expressions.html#StructuredOptimization.operator",
    "page": "Expressions",
    "title": "StructuredOptimization.operator",
    "category": "function",
    "text": "operator(ex::Expression)\n\nReturns the AbstractOperator of expression ex.\n\nExample\n\njulia> x = Variable(3)\nVariable(Float64, (3,))\n\njulia> ex = fft(x);\n\njulia> operator(ex)\nℱ  ℝ^3 -> ℂ^3\n\n\n\n\n"
},

{
    "location": "expressions.html#StructuredOptimization.displacement",
    "page": "Expressions",
    "title": "StructuredOptimization.displacement",
    "category": "function",
    "text": "displacement(ex::Expression)\n\nReturns the displacement of expression ex.\n\nExample\n\njulia> x = Variable(3)\nVariable(Float64, (3,))\n\njulia> ex = fft(x)+[1.+im*2.;0.;3.+im*4];\n\njulia> displacement(ex)\n3-element Array{Complex{Float64},1}:\n 1.0+2.0im\n 0.0+0.0im\n 3.0+4.0im\n\n\n\n\n"
},

{
    "location": "expressions.html#Utilities-2",
    "page": "Expressions",
    "title": "Utilities",
    "category": "section",
    "text": "It is possible to access the variables, mappings and displacement of an expression.  Notice that these commands work also for the Terms described in Functions and constraints.variables\noperator\ndisplacement"
},

{
    "location": "functions.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "functions.html#Functions-and-constraints-1",
    "page": "Functions",
    "title": "Functions and constraints",
    "category": "section",
    "text": "Once an expression is created it is possible to create the Terms defining the optimization problem. These can consists of either Smooth functions,  Nonsmooth functions, Inequality constraints  or Equality constraints."
},

{
    "location": "functions.html#StructuredOptimization.ls",
    "page": "Functions",
    "title": "StructuredOptimization.ls",
    "category": "function",
    "text": "ls(x::AbstractExpression)\n\nReturns the squared norm (least squares) of x:\n\nf (mathbfx) = frac12  mathbfx ^2\n\n(shorthand of 1/2*norm(x)^2).\n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.huberloss",
    "page": "Functions",
    "title": "StructuredOptimization.huberloss",
    "category": "function",
    "text": "huberloss(x::AbstractExpression, ρ=1.0)\n\nApplies the Huber loss function: \n\nf(mathbfx) = begincases\n  tfrac12 mathbfx ^2  textif  mathbfx  leq rho \n  rho ( mathbfx  - tfracrho2)  textotherwise\nendcases\n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.sqrhingeloss",
    "page": "Functions",
    "title": "StructuredOptimization.sqrhingeloss",
    "category": "function",
    "text": "sqrhingeloss(x::AbstractExpression, y::Array)\n\nApplies the squared Hinge loss function \n\nf( mathbfx ) = sum_i max0 1 - y_i x_i ^2\n\nwhere y is an array containing y_i.\n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.crossentropy",
    "page": "Functions",
    "title": "StructuredOptimization.crossentropy",
    "category": "function",
    "text": "crossentropy(x::AbstractExpression, y::Array)\n\nApplies the cross entropy loss function: \n\nf(mathbfx) = -1N sum_i^N y_i log (x_i)+(1-y_i) log (1-x_i)\n\nwhere y is an array of length N containing y_i having 0 leq y_i leq 1.\n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.logisticloss",
    "page": "Functions",
    "title": "StructuredOptimization.logisticloss",
    "category": "function",
    "text": "logbarrier(x::AbstractExpression, y::AbstractArray)\n\nApplies the logistic loss function: \n\nf(mathbfx) = sum_i log(1+ exp(-y_i x_i)) \n\nwhere y is an array containing y_i.\n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.logbarrier",
    "page": "Functions",
    "title": "StructuredOptimization.logbarrier",
    "category": "function",
    "text": "logbarrier(x::AbstractExpression)\n\nApplies the log barrier function: \n\nf(mathbfx) = -sum_i log( x_i )\n\n\n\n"
},

{
    "location": "functions.html#Smooth-functions-1",
    "page": "Functions",
    "title": "Smooth functions",
    "category": "section",
    "text": "ls\nhuberloss\nsqrhingeloss\ncrossentropy\nlogisticloss\nlogbarrier"
},

{
    "location": "functions.html#Base.LinAlg.norm",
    "page": "Functions",
    "title": "Base.LinAlg.norm",
    "category": "function",
    "text": "norm(x::AbstractExpression, p=2, [q,] [dim=1])\n\nReturns the norm of x. \n\nSupported norms:\n\np = 0 l_0-pseudo-norm\np = 1 l_1-norm\np = 2 l_2-norm\np = Inf l_infty-norm\np = * nuclear norm\np = 2, q = 1 l_21 mixed norm (aka Sum-of-l_2-norms) \n\nf(mathbfX) = sum_i  mathbfx_i \n\nwhere mathbfx_i is the i-th column if dim == 1 (or row if  dim == 2) of mathbfX.\n\n\n\n"
},

{
    "location": "functions.html#Base.maximum",
    "page": "Functions",
    "title": "Base.maximum",
    "category": "function",
    "text": "maximum(x::AbstractExpression)\n\nApplies the function: \n\nf(mathbfx) = max x_i  i = 1ldots n \n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.sumpositive",
    "page": "Functions",
    "title": "StructuredOptimization.sumpositive",
    "category": "function",
    "text": "sumpositive(x::AbstractExpression, ρ=1.0)\n\nApplies the function: \n\nf(mathbfx) = sum_i max x_i 0\n\n\n\n"
},

{
    "location": "functions.html#StructuredOptimization.hingeloss",
    "page": "Functions",
    "title": "StructuredOptimization.hingeloss",
    "category": "function",
    "text": "hingeloss(x::AbstractExpression, y::Array)\n\nApplies the Hinge loss function \n\nf( mathbfx ) = sum_i max0 1 - y_i x_i \n\nwhere y is an array containing y_i.\n\n\n\n"
},

{
    "location": "functions.html#Nonsmooth-functions-1",
    "page": "Functions",
    "title": "Nonsmooth functions",
    "category": "section",
    "text": "norm\nmaximum\nsumpositive\nhingeloss"
},

{
    "location": "functions.html#Base.:<=",
    "page": "Functions",
    "title": "Base.:<=",
    "category": "function",
    "text": "Inequalities constrains \n\nNorm Inequalities constraints\n\nnorm(x::AbstractExpression, 0) <= n::Integer \nmathrmnnz(mathbfx) leq n\nnorm(x::AbstractExpression, 1) <= r::Number \nsum_i  x_i  leq r\nnorm(x::AbstractExpression, 2) <= r::Number \n mathbfx  leq r\nnorm(x::AbstractExpression, Inf) <= r::Number \nmaxx_1 x_2 dots  leq r\n\nBox inequality constraints\n\nx::AbstractExpression <= u::Union{AbstractArray, Real}\nx_i leq u_i\nx::AbstractExpression >= l::Union{AbstractArray, Real}\nx_i geq l_i\nNotice that the notation x in [l,u] is also possible.\n\nRank inequality constraints\n\nrank(X::AbstractExpression) <= n::Integer\nmathrmrank(mathbfX) leq r \nNotice that the expression X must have a codomain with dimension equal to 2. \n\n\n\n"
},

{
    "location": "functions.html#Inequality-constraints-1",
    "page": "Functions",
    "title": "Inequality constraints",
    "category": "section",
    "text": "<="
},

{
    "location": "functions.html#Base.:==",
    "page": "Functions",
    "title": "Base.:==",
    "category": "function",
    "text": "Equalities constraints\n\nAffine space constraint\n\nex == b::Union{Real,AbstractArray} \nRequires expression to be affine.\nExample\njulia> A,b  = randn(10,5), randn(10);\n\njulia> x = Variable(5);\n\njulia> A*x == b\n\nNorm equality constraint\n\nnorm(x::AbstractExpression) == r::Number \n mathbfx  = r\n\nBinary constraint\n\nx::AbstractExpression == (l, u) \nmathbfx = mathbfl or mathbfx = mathbfu\n\n\n\n"
},

{
    "location": "functions.html#Equality-constraints-1",
    "page": "Functions",
    "title": "Equality constraints",
    "category": "section",
    "text": "=="
},

{
    "location": "functions.html#StructuredOptimization.smooth",
    "page": "Functions",
    "title": "StructuredOptimization.smooth",
    "category": "function",
    "text": "smooth(t::Term, gamma = 1.0)\n\nSmooths the nonsmooth term t using Moreau envelope:\n\nf^gamma(mathbfx) = min_mathbfz left f(mathbfz) + tfrac12gammamathbfz-mathbfx^2 right\n\nExample\n\njulia> x = Variable(4)\nVariable(Float64, (4,))\n\njulia> x = Variable(4);\n\njulia> t = smooth(norm(x,1))\n\n\n\n\n"
},

{
    "location": "functions.html#Smoothing-1",
    "page": "Functions",
    "title": "Smoothing",
    "category": "section",
    "text": "Sometimes the optimization problem might involve only non-smooth terms which do not lead to efficient proximal mappings. It is possible to smooth this terms by means of the Moreau envelope.smooth"
},

{
    "location": "functions.html#Base.conj",
    "page": "Functions",
    "title": "Base.conj",
    "category": "function",
    "text": "conj(t::Term)\n\nReturns the convex conjugate transform of t:\n\nf^*(mathbfx) = sup_mathbfy  langle mathbfy mathbfx rangle - f(mathbfy) \n\nExample\n\njulia> x = Variable(4)\nVariable(Float64, (4,))\n\njulia> x = Variable(4);\n\njulia> t = conj(norm(x,1))\n\n\n\n\n"
},

{
    "location": "functions.html#Duality-1",
    "page": "Functions",
    "title": "Duality",
    "category": "section",
    "text": "In some cases it is more convenient to solve the dual problem instead of the primal problem. It is possible to convert the primal problem into its dual form by means of the convex conjugate. See the Total Variation demo for an example of such procedure.conj"
},

{
    "location": "solvers.html#",
    "page": "Solvers",
    "title": "Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers.html#Solvers-1",
    "page": "Solvers",
    "title": "Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "solvers.html#StructuredOptimization.@minimize",
    "page": "Solvers",
    "title": "StructuredOptimization.@minimize",
    "category": "macro",
    "text": "@minimize cost [st ctr] [with slv_opt]\n\nMinimize a given problem with cost function cost, constraints ctr and solver options slv_opt. \n\nExample\n\njulia> A, b = randn(10,4), randn(10);\n\njulia> @minimize ls(A*x-b) + 0.5*norm(x);\n    it |      gamma |        fpr |        tau |        FBE |\n ------|------------|------------|------------|------------|\n     1 | 2.9152e-02 | 2.7656e+00 | 1.0000e+00 | 5.5181e+00 |\n     9 | 2.9152e-02 | 9.9682e-05 | 1.0000e+00 | 4.4086e+00 |\n\njulia> @minimize ls(A*x-b) st x >= 0.;\n    it |      gamma |        fpr |        tau |        FBE |\n ------|------------|------------|------------|------------|\n     1 | 5.8304e-02 | 1.0068e+00 | 1.0000e+00 | 6.6282e+00 |\n     3 | 5.8304e-02 | 9.5210e-16 | 1.0000e+00 | 6.5654e+00 |\n\njulia> it, slv = @minimize ls(A*x-b) st norm(x) == 2.0 with PG(maxit = 5);\n    it |      gamma |        fpr |\n ------|------------|------------|\n     1 | 6.1373e-02 | 2.2090e+01 |\n     5 | 3.0686e-02 | 5.5190e-01 |\n\n\nReturns as output a tuple containing the number of iterations and the constructed solver.\n\n\n\n"
},

{
    "location": "solvers.html#Minimizing-a-problem-1",
    "page": "Solvers",
    "title": "Minimizing a problem",
    "category": "section",
    "text": "@minimizenote: Problem warm-starting\nBy default warm-starting is always enabled. For example, if two problems that utilize the same variables are solved consecutively,  the second one will be automatically warm-started by the solution of the first one.That is because the variables are always linked to their respective data vectors. If one wants to avoid this, the optimization variables needs to be manually re-initialized  before solving the second problem e.g. to a vector of zeros: ~x .= 0.0."
},

{
    "location": "solvers.html#StructuredOptimization.PG",
    "page": "Solvers",
    "title": "StructuredOptimization.PG",
    "category": "type",
    "text": "PG(;kwargs...)\n\nCreates an object PG containing the options of the Proximal Gradient solvers:\n\ngamma, stepsize (default: unspecified, determined automatically)\nmaxit, maximum number of iteration (default: 10000)\ntol, halting tolerance on the fixed-point residual (default: 1e-4)\nadaptive, adaptively adjust gamma (default: false if gamma is provided)\nfast, enables accelerated method (default: false)\nverbose, verbosity level (default: 1)\nverbose_freq, verbosity frequency for verbose = 1 (default: 100)\n\n\n\n"
},

{
    "location": "solvers.html#StructuredOptimization.FPG",
    "page": "Solvers",
    "title": "StructuredOptimization.FPG",
    "category": "function",
    "text": "FPG(;kwargs...)\n\nSame as PG, creates the options of the Fast Proximal Gradient solver. \n\n\n\n"
},

{
    "location": "solvers.html#StructuredOptimization.ZeroFPR",
    "page": "Solvers",
    "title": "StructuredOptimization.ZeroFPR",
    "category": "type",
    "text": "ZeroFPR(;kwargs...)\n\nCreates an object ZeroFPR containing the options of the ZeroFPR solver:\n\ngamma, stepsize (default: unspecified, determined automatically)\nmaxit, maximum number of iteration (default: 10000)\ntol, halting tolerance on the fixed-point residual (default: 1e-4)\nadaptive, adaptively adjust gamma (default: false if gamma is provided)\nfast, enables accelerated method (default: false)\nverbose, verbosity level (default: 1)\nverbose_freq, verbosity frequency for verbose = 1 (default: 100)\nmemory, memory of the LBFGS operator (default: 10 )\n\n\n\n"
},

{
    "location": "solvers.html#StructuredOptimization.PANOC",
    "page": "Solvers",
    "title": "StructuredOptimization.PANOC",
    "category": "type",
    "text": "ZeroFPR(;kwargs...)\n\nCreates an object PANOC containing the options of the PANOC solver:\n\ngamma, stepsize (default: unspecified, determined automatically)\nmaxit, maximum number of iteration (default: 10000)\ntol, halting tolerance on the fixed-point residual (default: 1e-4)\nadaptive, adaptively adjust gamma (default: false if gamma is provided)\nfast, enables accelerated method (default: false)\nverbose, verbosity level (default: 1)\nverbose_freq, verbosity frequency for verbose = 1 (default: 100)\nmemory, memory of the LBFGS operator (default: 10 )\n\n\n\n"
},

{
    "location": "solvers.html#Specifying-solver-and-options-1",
    "page": "Solvers",
    "title": "Specifying solver and options",
    "category": "section",
    "text": "As shown above it is possible to choose the type of algorithm and specify its options by creating a Solver object.Currently, the following algorithms are supported:Proximal Gradient (PG) [1], [2]\nFast Proximal Gradient (FPG) [1], [2]\nZeroFPR [3]\nPANOC [4]PG\nFPG\nZeroFPR\nPANOC"
},

{
    "location": "solvers.html#StructuredOptimization.problem",
    "page": "Solvers",
    "title": "StructuredOptimization.problem",
    "category": "function",
    "text": "problems(terms...)\n\nConstructs a problem.\n\nExample\n\n\njulia> x = Variable(4)\nVariable(Float64, (4,))\n\njulia> A, b = randn(10,4), randn(10);\n\njulia> p = problem(ls(A*x-b), norm(x) <= 1)\n\n\n\n\n"
},

{
    "location": "solvers.html#StructuredOptimization.solve",
    "page": "Solvers",
    "title": "StructuredOptimization.solve",
    "category": "function",
    "text": "solve(terms::Tuple, solver_opt::ForwardBackwardSolver)\n\nTakes as input a tuple containing the terms defining the problem and the solver options.\n\nSolves the problem returning a tuple containing the iterations taken and the build solver.\n\nExample\n\n\njulia> x = Variable(4)\nVariable(Float64, (4,))\n\njulia> A, b = randn(10,4), randn(10);\n\njulia> solve(p,PG());\n    it |      gamma |        fpr |\n ------|------------|------------|\n     1 | 7.6375e-02 | 1.8690e+00 |\n    12 | 7.6375e-02 | 9.7599e-05 |\n\n\n\n\n"
},

{
    "location": "solvers.html#StructuredOptimization.build",
    "page": "Solvers",
    "title": "StructuredOptimization.build",
    "category": "function",
    "text": "build(terms::Tuple, solver_opt::ForwardBackwardSolver)\n\nTakes as input a tuple containing the terms defining the problem and the solver options.\n\nReturns a tuple containing the optimization variables and the built solver.\n\nExample\n\njulia> x = Variable(4)\nVariable(Float64, (4,))\n\njulia> A, b = randn(10,4), randn(10);\n\njulia> p = problem( ls(A*x - b ) , norm(x) <= 1 );\n\njulia> build(p, PG());\n\n\n\n\n"
},

{
    "location": "solvers.html#StructuredOptimization.solve!",
    "page": "Solvers",
    "title": "StructuredOptimization.solve!",
    "category": "function",
    "text": "solve!( x_solver )\n\nTakes as input a tuple containing the optimization variables and the built solver.\n\nSolves the problem returning a tuple containing the iterations taken and the build solver.\n\nExample\n\njulia> x = Variable(4)\nVariable(Float64, (4,))\n\njulia> A, b = randn(10,4), randn(10);\n\njulia> p = problem( ls(A*x - b ) , norm(x) <= 1 );\n\njulia> x_solver = build(p, PG(verbose = 0));\n\njulia> solve!(x_solver);\n\n\n\n\n"
},

{
    "location": "solvers.html#Build-and-solve-1",
    "page": "Solvers",
    "title": "Build and solve",
    "category": "section",
    "text": "The macro @minimize automatically parse and solve the problem. An alternative syntax is given by the function problem and solve.problem\nsolveIt is important to stress out that the Solver objects created using the functions above (PG, FPG, etc.) specify only the type of algorithm to be used together with its options. The actual solver  (namely the one of ProximalAlgorithms.jl)  is constructed altogether with the problem formulation. The problem parsing procedure can be separated from the solver application using the functions build and solve!.build\nsolve!"
},

{
    "location": "solvers.html#References-1",
    "page": "Solvers",
    "title": "References",
    "category": "section",
    "text": "[1] Tseng, On Accelerated Proximal Gradient Methods for Convex-Concave Optimization (2008).[2] Beck, Teboulle, A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems, SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183-202 (2009).[3] Themelis, Stella, Patrinos, Forward-backward envelope for the sum of two nonconvex functions: Further properties and nonmonotone line-search algorithms, arXiv:1606.06256 (2016).[4] Stella, Themelis, Sopasakis, Patrinos, A simple and efficient algorithm for nonlinear model predictive control, 56th IEEE Conference on Decision and Control (2017)."
},

{
    "location": "demos.html#",
    "page": "Demos",
    "title": "Demos",
    "category": "page",
    "text": ""
},

{
    "location": "demos.html#Demos-1",
    "page": "Demos",
    "title": "Demos",
    "category": "section",
    "text": "Sparse deconvolution\nLine Spectra Estimation\nDeep neural network classifier\nVideo background removal\nTotal variation denoising\nAudio declippingClipped audio sample (Warning there are severe distortions and you might want to turn down your volume before playing)<audio src=\"assets/clipped.wav\" controls preload></audio>De-clipped audio sample <audio src=\"assets/declipped.wav\" controls preload></audio>"
},

]}
