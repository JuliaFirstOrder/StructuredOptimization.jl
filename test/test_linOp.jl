##test linear operators

stuff = [
#	 Dict("Operator" => (zeros,),
#              "params" => ((),),
#	      "args"   => ( randn(4), randn(4) )
#	     ),
#	 Dict("Operator" => (zeros,),
#              "params" => ((2,),),
#	      "args"   => ( randn(4), randn(2) )
#	     ),
	 Dict(
       "Operator" => (Eye,),
       "params"   => ((Float64, (4,4),),),
       "args"     => ( randn(4,4), randn(4,4) )
	     ),
#	 Dict("Operator" => (diagop,),
#              "params" => ((randn(2,2)+im*randn(2,2),),),
#	      "args"   => ( randn(2,2)+im*randn(2,2), randn(2,2)+im*randn(2,2) )
#	      ),
#	 Dict("Operator" => (diagop,),
#              "params" => ((2,),),
#	      "args"   => ( randn(2,2), randn(2,2) )
#	      ),
#	 Dict("Operator" => (diagop,),
#              "params" => ((2+im*3,),),
#	      "args"   => ( randn(2,2), (2+im*3)*randn(2,2) )
#	      ),
#	 Dict("Operator" => (*,),
#              "params" => ((2+im*3),),
#	      "args"   => ( randn(2,2), (2+im*3)*randn(2,2) )
#	      ),
#	 Dict("Operator" => (.*,),
#              "params" => ((randn(2,2)),),
#	      "args"   => ( randn(2,2), randn(2,2) )
#	      ),
#	 Dict("Operator" => (*,),
#              "params" => ((randn(4,4)),),
#	      "args"   => ( randn(4), randn(4) )
#	      ),
#	 Dict("Operator" => (*,),
#              "params" => ((randn(6,2)+im*randn(6,2)),),
#	      "args"   => ( randn(2)+im*randn(2), randn(6)+im*randn(6) )
#	      ),
#	 Dict("Operator" => (*,),
#              "params" => ((randn(2,6)+im*randn(2,6)),),
#	      "args"   => ( randn(6)+im*randn(6), randn(2)+im*randn(2) )
#	      ),
#	 Dict("Operator" => (getindex,),
#              "params" => (([1:3]),),
#	      "args"   => ( randn(6)+im*randn(6), randn(3)+im*randn(3) )
#	      ),
#	 Dict("Operator" => (getindex,),
#              "params" => (([1:3,:]),),
#	      "args"   => ( randn(4,2)+im*randn(4,2), randn(3,2)+im*randn(3,2) )
#	      ),
#	 Dict("Operator" => (getindex,),
#              "params" => (([1:3]),),
#	      "args"   => ( randn(4,2)+im*randn(4,2), randn(3)+im*randn(3) )
#	      ),
#	 Dict("Operator" => (fft,),
#              "params" => ((), ),
#	      "args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
#	      ),
#	 Dict("Operator" => (fft,),
#              "params" => ((),),
#	      "args"   => ( randn(64,64), fft(randn(64,64)) )
#	      ),
#	 Dict("Operator" => (ifft,),
#              "params" => ((), ),
#	      "args"   => ( randn(64,64), fft(randn(64,64)) )
#	      ),
#	 Dict("Operator" => (ifft,),
#              "params" => ((), ),
#	      "args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
#	      ),
#	 Dict("Operator" => (dct,),
#	      "params" => ((), ),
#	      "args"   => ( randn(64,64), randn(64,64) )
#	      ),
#	 Dict("Operator" => (idct,),
#              "params" => ((), ),
#	      "args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
#	      ),
#
#	 Dict("Operator" => (reshape,),
#              "params" => ((10,10), ),
#	      "args"   => ( randn(100), randn(10,10) )
#	      ),
#	 Dict("Operator" => (conv,),
#              "params" => (([randn(100)]),),
#	      "args"   => ( randn(400), randn(400+100-1) )
#	     ),
#	 Dict("Operator" => (xcorr,),
#              "params" => (([ones(100)]),),
#	      "args"   => ( ones(30), ones(2*100-1) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((),),
#	      "args"   => ( randn(4), randn(4) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((),),
#	      "args"   => ( randn(4), randn(4) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((1),),
#	      "args"   => ( randn(4,6), randn(4,6) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((2),),
#	      "args"   => ( randn(4,6), randn(4,6) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((1),),
#	      "args"   => ( randn(4,6,3), randn(4,6,3) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((2),),
#	      "args"   => ( randn(4,6,3), randn(4,6,3) )
#	     ),
#	 Dict("Operator" => (finiteDiff,),
#              "params" => ((3),),
#	      "args"   => ( randn(4,6,3), randn(4,6,3) )
#	     ),
#	 Dict("Operator" => (tv,),
#              "params" => ((),),
#	      "args"   => ( randn(4,6), randn(4*6,2) )
#	     ),
#	 Dict("Operator" => (tv,),
#              "params" => ((),),
#	      "args"   => ( randn(4,6,5), randn(4*6*5,3) )
#	     ),
#	 #####some tests on nested operators###
#	 Dict("Operator" => (*, dct),
#              "params" => ((randn(32,64)) ,() ),
#	      "args"   => ( randn(64), randn(32) )
#	      ),
#	 Dict("Operator" => (ifft, reshape, dct),
#              "params" => ((),(10,10),() ),
#	      "args"   => ( randn(100), dct(reshape(ifft(randn(100)),10,10)) )
#	      ),
#	 Dict("Operator" => (dct, getindex, reshape, dct),
#              "params" => ((),([1:100]),(10,10),() ),
#	      "args"   => ( randn(120), randn(10,10) )
#	      ),
#	 Dict("Operator" => (ifft, getindex),
#              "params" => ((),([1:20]) ),
#	      "args"   => ( randn(200)+im*randn(200), randn(20)+im*randn(20) )
#	      ),
#	 Dict("Operator" => (ifft, getindex),
#              "params" => ((),([1:20]) ),
#	      "args"   => ( randn(10,10)+im*randn(10,10), randn(20)+im*randn(20) )
#	      ),
#	 Dict("Operator" => (fft, getindex),
#              "params" => ((),([1:5]) ),
#	      "args"   => ( randn(12)+im*randn(12), randn(5)+im*randn(5) )
#	      ),
#	 Dict("Operator" => (dct, getindex),
#              "params" => ((),([1:2,:,2:5]) ),
#	      "args"   => ( randn(5,5,5)+im*randn(5,5,5), randn(2,5,4)+im*randn(2,5,4) )
#	      ),
	 ]


for i in eachindex(stuff)

	x,y = deepcopy(stuff[i]["args"])

	params = stuff[i]["params"][1]
	Op     = stuff[i]["Operator"][1]
	A = Op(params...)
	# Aminus = -operator(A)
	# @test vecnorm(operator(A)*x + Aminus*x) <= 1e-8 #test sign
	test1,test2 = RegLS.test_FwAdj(A, x, y)
	@test test1 < 1e-8
	@test test2 < 1e-8
	test3 = RegLS.test_Op(A, x, y)
	@test test3 < 1e-8

end


##### test sum of linear operators

#test SumSameVar
#
#x1 = randn(3)
#y1 = randn(3)
#b1 = randn(3)
#M = randn(3,3)
#X1 = Variable(x1)
#x = x1
#y = y1
#
#A = -M*X1
#@test norm(A(x1)-(-M*x1)) <= 1e-8
#A = b1-M*X1
#@test norm(A(x1)-(b1-M*x1)) <= 1e-8
#A = X1+b1
#@test norm(A(x1)-(b1+x1)) <= 1e-8
#A = X1-b1
#@test norm(A(x1)-(-b1+x1)) <= 1e-8
#A = -X1+b1
#@test norm(A(x1)-(b1-x1)) <= 1e-8
#
#A = -M*X1
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
#x1 = randn(3)
#X1 = Variable(x1)
#y = randn(3)
#M = randn(3,3)
#x = x1
#b1 = randn(3)
#
#A = -3*X1-b1-dct(X1)
#
#@test norm(A(x)-(-3*x-b1-dct(x))) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
#A = -3*X1-(b1-(dct(X1)+M*X1))
#show(A)
#
#@test norm(A(x)-(-3*x-(b1-(dct(x)+M*x)))) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
#A = -3*X1-b1-dct(X1)-M*X1
#
#@test norm(A(x)-(-3*x-b1-dct(x)-M*x)) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
###test HCAT
#
#x1,x2 = randn(3,3), randn(3,3)
#X1,X2 = Variable(x1), Variable(x2)
#y = randn(3,3)
#b = randn(3,3)
#x = [x1,x2]
#
#A = 3.4*X1-2.0*dct(X2)+b
#A2 = [DiagOp(3.4,3,3)*Eye(3,3) -2.0*DCT(X2.x)]
#@test norm(A(x)-(3.4*x1-2*dct(x2)+b)) < 1e-8
#@test norm(A(x)-(A2*x+b)) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
###test HCAT .*
#
#x1,x2 = randn(3,3), randn(3*3)
#X1,X2 = Variable(x1), Variable(x2)
#y = randn(3,3)
#x = [x1,x2]
#
#B = [Eye(3,3) Reshape(X2.x,3,3)*DCT(X2.x)]
#A = [DiagOp(3.4,3,3) DiagOp(-2.0,3,3) ]
#C = Affine([X1,X2],A.*B)
#
#@test norm(C(x)-(3.4*x1-2.0*reshape(dct(x2),3,3))) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(C, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(C, x, y)
#@test test3 < 1e-8
#
#A.*[x1,x1]
#
##test HCAT merge
#
#x1,x2,x3 = randn(3),randn(3),randn(3)
#X1,X2,X3 = Variable(x1), Variable(x2), Variable(x3)
#M = randn(3,3)
#y = randn(3)
#b = randn(3)
#x = [x1,x2,x3]
#
#A = M*X1-X1+5.3*(M*X2)-b+eye(X1)+X3
#A2 = [MatrixOp(M) 5.3*MatrixOp(M) Eye(3)]
#@test norm(A(x)-(M*x1-x1+5.3*M*x2+x1+x3-b)) < 1e-8
#@test norm((A2*x)-(M*x1-x1+5.3*M*x2+x1+x3)) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
#A1 = M*X1-X2
#B1 = M*X3+4.4*X2
#A = A1-B1
#
#@test norm(A(x)-(M*x1-x2-(M*x3+4.4*x2))) < 1e-8
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
##test VCAT
#
#x1 = randn(3,3)
#X1 = Variable(x1)
#y1,y2 = randn(3,3),randn(3)
#Y1,Y2 = Variable(y1),Variable(y2)
#x = x1
#y = [y1,y2]
#
#A = [operator(3.4*X1); operator(-2.0*dct(X1)[1:3])]
#show(A)
#
#@test norm((A*x)[1]-3.4*x1) < 1e-8
#@test norm((A*x)[2]-(-2*dct(x1)[1:3])) < 1e-8
#
#A = Affine([Y1,Y2],A')
#
#test1,test2 = RegLS.test_FwAdj(A, y, x)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, y, x)
#@test test3 < 1e-8
#
##test VCAT .*
#
#x1 = randn(3,3)
#X1 = Variable(x1)
#y1,y2 = randn(3,3),randn(6)
#x = x1
#y = [y1,y2]
#
#A = [operator(5*X1);    operator(reshape(5*X1[1:2,:],6)) ]
#B = [operator(dct(X1)); operator(4*dct(X1))  ]
#C = A.*B
#
#@test norm((C*x)[1]-5.*dct(x1))<=1e-8
#@test norm((C*x)[2]-reshape((5*(4*dct(x1))[1:2,:]),6) )<=1e-8
#
#C = Affine(X1,C')
#
#test1,test2 = RegLS.test_FwAdj(C, y, x)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(C, y, x)
#@test test3 < 1e-8
#
#
###test sorting Affine
#
#x,y = Variable(10), Variable(2,5)
#X = [randn(10), randn(2,5)]
#Y = randn(10)
#A = eye(x)-reshape(y,10)
#
#show(issorted(A))
#if issorted(A) == false #this if is in case object_id are actually sorted already!
#
#	show(operator(A).A)
#	y1 = operator(A)*X
#	x1 = adjoint(A)*Y
#	p = sortperm(A)
#	A2 = sort(A)
#	X2 = X[p]
#	y2 = operator(A2)*X2
#	x2 = adjoint(A2)*Y
#
#	show(operator(A2).A)
#	@test issorted(A2) == true
#	@test norm(y1-y2)<=1e-8
#	@test norm(x1[p]-x2)<=1e-8
#else
#	A = -reshape(y,10)+eye(x)
#	X = [X[2], X[1]]
#
#	show(operator(A).A)
#	y1 = operator(A)*X
#	x1 = adjoint(A)*Y
#	p = sortperm(A)
#	A2 = sort(A)
#	X2 = X[p]
#	y2 = operator(A2)*X2
#	x2 = adjoint(A2)*Y
#
#	show(operator(A2).A)
#	@test issorted(A2) == true
#	@test norm(y1-y2)<=1e-8
#	@test norm(x1[p]-x2)<=1e-8
#
#end
#
###test sort_and_expand
#
#x,y,z = Variable(10), Variable(2,5), Variable(10)
#X = [randn(10), randn(2,5), randn(10)]
#Y = randn(10)
#A = reshape(y,10)
#show(A)
#
#B = RegLS.sort_and_expand([x,y,z],A)
#
#test1,test2 = RegLS.test_FwAdj(B, X, Y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(B, X, Y)
#@test test3 < 1e-8
#
#A = reshape(y,10)+2.5*x
#show(A)
#
#B = RegLS.sort_and_expand([x,y,z],A)
#
#test1,test2 = RegLS.test_FwAdj(B, X, Y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(B, X, Y)
#@test test3 < 1e-8
#
#A = reshape(y,10)+2.5*x+randn(10,10)*z
#show(A)
#
#B = RegLS.sort_and_expand([x,y,z],A)
#
#test1,test2 = RegLS.test_FwAdj(B, X, Y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(B, X, Y)
#@test test3 < 1e-8
#
##testing FiniteDiff
#
#b = repmat(collect(linspace(0,1,5)),1,4)
#x = Variable(b)
#
#A = FiniteDiff(~x)
#
#y = A*b
#@test vecnorm(y-0.25*ones(size(b))) <1e-8
#
#A = FiniteDiff(~x, 2)
#y = A*b
#@test vecnorm(y) <1e-8
#
#A = TV(~x)
#y = A*b
#
#@test vecnorm(y[:,1]-0.25*ones(size(b))[:]) <1e-8
#@test vecnorm(y[:,2]) <1e-8
#
#b = randn(size(b))
#y1 = TV(~x)*b
#yx = FiniteDiff(~x)*b
#yy = FiniteDiff(~x,2)*b
#
#@test norm(y1[:,1]-yx[:]) <=1e-8
#@test norm(y1[:,2]-yy[:]) <=1e-8
#
#b = randn(40,50)
#x = Variable(b)
#
#x1 = FiniteDiff(~x)'*b
#x2 = TV(~x)'*[b[:] 0*b[:]]
#
#@test vecnorm(x1-x2) <= 1e-8
#
#x1 = FiniteDiff(~x,2)'*b
#x2 = TV(~x)'*[0*b[:] b[:]]
#
#@test vecnorm(x1-x2) <= 1e-8
#
#b = randn(40,50,20)
#x = Variable(b)
#
#x1 = FiniteDiff(~x)'*b
#x2 = TV(~x)'*[b[:] 0*b[:] 0*b[:]]
#
#@test vecnorm(x1-x2) <= 1e-8
