##test linear operators

stuff = [
### testing constructors ###
	 Dict(
       "Operator" => (Zeros,),
       "params"   => ((4,4,),),
       "args"     => ( randn(4,4), randn(4,4) ),
       "in_out"   => ( ones(4,4), zeros(4,4) )
	     ),
	 Dict(
       "Operator" => (Zeros,),
       "params"   => ((randn(4,4),),),
       "args"     => ( randn(4,4), randn(4,4) ),
       "in_out"   => ( ones(4,4), zeros(4,4) )
	     ),
	 Dict(
       "Operator" => (Zeros,),
       "params"   => ((MatrixOp(randn(4,3)),),),
       "args"     => ( randn(3), randn(4) ),
	     ),
	 Dict(
       "Operator" => (Eye,),
       "params"   => ((4,4,),),
       "args"     => ( randn(4,4), randn(4,4) ),
       "in_out"   => ( ones(4,4), ones(4,4) )
	     ),
	 Dict(
       "Operator" => (Eye,),
       "params"   => ((Complex{Float64}, (4,4),),),
       "args"     => ( randn(4,4)+im, randn(4,4)+im ),
       "in_out"   => ( ones(4,4)+im*ones(4,4), ones(4,4)+im*ones(4,4) )
	     ),
	 Dict(
       "Operator" => (GetIndex,),
       "params"   => (( zeros(4,4), (1:2,:)  ,),),
       "args"     => ( randn(4,4), randn(2,4) ),
       "in_out"   => ( ones(4,4), ones(2,4) )
	     ),
	 Dict(
       "Operator" => (GetIndex,),
       "params"   => (( (4,4), (1:2,:)  ,),),
       "args"     => ( randn(4,4), randn(2,4) ),
       "in_out"   => ( ones(4,4), ones(2,4) )
	     ),
	 Dict(
       "Operator" => (GetIndex,),
       "params"   => ((Complex{Float64}, (4,4), (1:2,:)  ,),),
       "args"     => ( randn(4,4)+im, randn(2,4)+im ),
       "in_out"   => ( ones(4,4)+im, ones(2,4)+im )
	     ),
	 Dict(
       "Operator" => (Eye,),
       "params"   => ((randn(4,4),),),
       "args"     => ( randn(4,4), randn(4,4) ),
	     ),
	 Dict(
       "Operator" => (DiagOp,),
       "params"   => ((randn(srand(1),4),),),
       "wrg_pr"   => ((randn(srand(1),4), randn(srand(1),4)+im,),),
       "args"     => ( randn(4), randn(4) ),
       "in_out"   => (ones(4),randn(srand(1),4).*ones(4))
	     ),
	 Dict(
       "Operator" => (DiagOp,),
       "params"   => ((randn(4,4)*im, randn(srand(1),4,4)+im*randn(srand(2),4,4),),),
       "wrg_pr"   => ((randn(4,4), randn(srand(1),4,4)+im*randn(srand(2),4,4),),),
       "args"     => ( randn(4,4)+randn(4,4)*im, randn(4,4)+randn(4,4)*im ),
	     ),
	 Dict(
       "Operator" => (MatrixOp,),
       "params"   => ((randn(srand(1),4,6),),),
       "args"     => ( randn(6), randn(4) ),
       "in_out"   => (ones(6),randn(srand(1),4,6)*ones(6))
	     ),
	 Dict(
       "Operator" => (MatrixOp,),
       "params"   => ((randn(4,6)+im*randn(4,6),),),
       "args"     => ( randn(6)+im, randn(4)+im )
	     ),
	 Dict(
       "Operator" => (MatrixOp,),
       "params"   => ((randn(4,6)+im*randn(4,6),),),
       "args"     => ( randn(6)+im, randn(4)+im )
	     ),
	 Dict(
       "Operator" => (MatrixOp,),
       "params"   => ((randn(4,3),3,2),),
       "args"     => ( randn(3,2), randn(4,2) )
	     ),
	 Dict(
       "Operator" => (MatrixOp,),
       "params"   => ((Complex{Float64}, randn(4,3), 3,2 ),),
       "args"     => ( randn(3,2)+im, randn(4,2)+im )
	     ),
	 Dict(
       "Operator" => (MatrixOp,),
       "params"   => ((randn(3,2)+im, randn(4,3)),),
       "args"     => ( randn(3,2)+im, randn(4,2)+im )
	     ),
	 Dict(
       "Operator" => (DFT,),
       "params"   => ((randn(3,2)+im,),),
       "args"     => ( randn(3,2)+im, randn(3,2)+im ),
       "in_out"   => ( ones(3,2)+2*im*ones(3,2), fft(ones(3,2)+2*im*ones(3,2)) )
	     ),
	 Dict(
       "Operator" => (DFT,),
       "params"   => ((10,),),
       "args"     => ( randn(10), fft(randn(10)) ),
	     ),
	 Dict(
       "Operator" => (DFT,),
       "params"   => ((Float64,10,),),
       "args"     => ( randn(10), fft(randn(10)) ),
	     ),
	 Dict(
       "Operator" => (IDFT,),
       "params"   => ((randn(3,2)+im,),),
       "args"     => ( randn(3,2)+im, randn(3,2)+im ),
       "in_out"   => ( ones(3,2)+2*im*ones(3,2), ifft(ones(3,2)+2*im*ones(3,2)) )
	     ),
	 Dict(
       "Operator" => (IDFT,),
       "params"   => ((10,),),
       "args"     => ( randn(10), fft(randn(10)) ),
	     ),
	 Dict(
       "Operator" => (IDFT,),
       "params"   => ((Float64,10),),
       "args"     => ( randn(10), fft(randn(10)) ),
	     ),
	 Dict(
       "Operator" => (DCT,),
       "params"   => ((randn(3,2)+im,),),
       "args"     => ( randn(3,2)+im, randn(3,2)+im ),
       "in_out"   => ( ones(3,2)+2*im*ones(3,2), dct(ones(3,2)+2*im*ones(3,2)) )
	     ),
	 Dict(
       "Operator" => (IDCT,),
       "params"   => ((randn(3,2),),),
       "args"     => ( randn(3,2), randn(3,2) ),
       "in_out"   => ( ones(3,2), idct(ones(3,2)) )
	     ),
	 Dict(
       "Operator" => (IDCT,),
       "params"   => ((3,2,),),
       "args"     => ( randn(3,2), randn(3,2) ),
       "in_out"   => ( ones(3,2), idct(ones(3,2)) )
	     ),
	 Dict(
       "Operator" => (Conv,),
       "params"   => ((10,randn(srand(1),20),),),
       "args"     => ( randn(10), randn(29) ),
       "in_out"   => ( randn(srand(3),10), conv(randn(srand(1),20),randn(srand(3),10)) )
	     ),
	 Dict(
       "Operator" => (Conv,),
       "params"   => ((randn(10),randn(srand(1),20),),),
       "args"     => ( randn(10), randn(29) ),
	     ),
	 Dict(
       "Operator" => (Xcorr,),
       "params"   => ((randn(30),randn(100),),),
       "args"     => ( randn(30), randn(199) ),
	     ),
	 Dict(
       "Operator" => (Xcorr,),
       "params"   => ((30,randn(100),),),
       "args"     => ( randn(30), randn(199) ),
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((randn(10),),),
       "args"     => ( randn(10), randn(10) ),
       "in_out"   => ( collect(linspace(0,1,10)), 1/9*ones(10) )
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => (((10,),),),
       "args"     => ( randn(10), randn(10) ),
       "in_out"   => ( collect(linspace(0,1,10)), 1/9*ones(10) )
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((Complex{Float64},(10,),),),
       "args"     => ( randn(10)+im, randn(10)+im ),
       "in_out"   => ( collect(linspace(0,1,10))+im*collect(linspace(0,1,10)), 
		      1/9*(ones(10)+im*ones(10)) )
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((randn(10,5),),),
       "args"     => ( randn(10,5), randn(10,5) ),
       "in_out"   => ( repmat(collect(linspace(0,1,10)),1,5), 1/9*ones(10,5) )
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((randn(10,5),2,),),
       "args"     => ( randn(10,5), randn(10,5) ),
       "in_out"   => ( repmat(collect(linspace(0,1,10)),1,5), zeros(10,5) )
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((randn(10,5,6),1,),),
       "args"     => ( randn(10,5,6), randn(10,5,6) ),
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((randn(10,5,6),2,),),
       "args"     => ( randn(10,5,6), randn(10,5,6) ),
	     ),
	 Dict(
       "Operator" => (FiniteDiff,),
       "params"   => ((randn(10,5,6),3,),),
       "args"     => ( randn(10,5,6), randn(10,5,6) ),
	     ),
	 Dict(
       "Operator" => (Variation,),
       "params"   => ((randn(10,5),),),
       "args"     => ( randn(10,5), randn(10*5,2) ),
       "in_out"   => ( repmat(collect(linspace(0,1,10)),1,5), [1/9*ones(10*5) zeros(10*5)] )
	     ),
	 Dict(
       "Operator" => (Variation,),
       "params"   => ((10,5,),),
       "args"     => ( randn(10,5), randn(10*5,2) ),
	     ),
	 Dict(
       "Operator" => (Variation,),
       "params"   => ((randn(10,5,4),),),
       "args"     => ( randn(10,5,4), randn(10*5*4,3) ),
	     ),
### testing Reshape ####
	 Dict(
       "Operator" => (Reshape,),
       "params"   => ((MatrixOp(2*ones(4,8)),2,2),),
       "args"     => ( randn(8), randn(2,2) ),
       "in_out"   => ( ones(8), 16*ones(2,2) )
	     ),
	 Dict(
       "Operator" => (Reshape,),
       "params"   => ((Variation(ones(5,4)),5*4*2),),
       "args"     => ( randn(5,4), randn(5*4*2) ),
	     ),
#### testing Scale ####
	 Dict(
       "Operator" => ((*),),
       "params"   => ((2, MatrixOp(randn(srand(1),4,8))),),
       "wrg_pr"   => ((2+im, MatrixOp(randn(4,8))),),
       "args"     => ( randn(8), randn(4) ),
       "in_out"   => ( randn(srand(3),8), 2*randn(srand(1),4,8)*randn(srand(3),8) )
	     ),
	 Dict(
       "Operator" => ((*),),
       "params"   => ((1+pi*im, DFT(randn(4)+im )),),
       "args"     => ( randn(4)+im, randn(4)+im ),
       "in_out"   => ( randn(srand(4),4)+im*randn(srand(5),4), 
		      (1+pi*im)*fft( randn(srand(4),4)+im*randn(srand(5),4)) )
	     ),
	 Dict(
       "Operator" => ((*),),
       "params"   => ((5, DFT(randn(4)+im )),),
       "args"     => ( randn(4)+im, randn(4)+im ),
       "in_out"   => ( randn(srand(4),4)+im*randn(srand(5),4), 
		      5*fft( randn(srand(4),4)+im*randn(srand(5),4)) )
	     ),
	 Dict(
       "Operator" => ((-),),
       "params"   => ((MatrixOp(randn(4,8)),),),
       "args"     => ( randn(8), randn(4) ),
	     ),
### testing Sum ####
	 Dict(
       "Operator" => ((+),),
       "params"   => ((MatrixOp(randn(srand(1),4,8)), MatrixOp(randn(srand(2),4,8))),),
       "wrg_pr"   => ((MatrixOp(randn(4,7)), MatrixOp(randn(4,8))),),
       "args"     => ( randn(8), randn(4) ),
       "in_out"   => ( randn(srand(3),8), 
		     (randn(srand(1),4,8)+randn(srand(2),4,8))*randn(srand(3),8) )
	     ),
	 Dict(
       "Operator" => ((-),),
       "params"   => ((2*DFT(randn(4)), DFT(randn(4)),),),
       "wrg_pr"   => ((DFT(randn(4)), MatrixOp(randn(4,4)),),),
       "args"     => ( randn(4), fft(randn(4)) ),
       "in_out"   => ( randn(srand(3),4), fft(randn(srand(3),4)))
	     ),
	 Dict(
       "Operator" => ((-),),
       "params"   => ((DFT(randn(4))+DFT(randn(4)),DFT(randn(4))),),
       "args"     => ( randn(4), fft(randn(4)) ),
       "in_out"   => ( randn(srand(3),4), fft(randn(srand(3),4)))
	     ),
	 Dict(
       "Operator" => ((-),),
       "params"   => ((DFT(randn(4)),DFT(randn(4))-DFT(randn(4)) ),),
       "args"     => ( randn(4), fft(randn(4)) ),
       "in_out"   => ( randn(srand(3),4), fft(randn(srand(3),4)))
	     ),
### testing Compose ####
	 Dict(
       "Operator" => ((*),),
       "params"   => (( MatrixOp(randn(srand(3),4,3)) , DCT((3,))                     ,),),
       "wrg_pr"   => (( DCT((3,))                     , MatrixOp(randn(srand(3),4,3)) ,),),
       "args"     => ( randn(3), randn(4) ),
       "in_out"   => ( randn(srand(3),3), randn(srand(3),4,3)*dct( randn(srand(3),3) )  )
	     ),
	 Dict(
       "Operator" => ((*),),
       "params"   => (( DFT((3,))                     , MatrixOp(randn(srand(3),3,4)) ,),),
       "wrg_pr"   => (( DFT(Complex{Float64},(3,))    , MatrixOp(randn(srand(3),3,4)) ,),),
       "args"     => ( randn(4), fft(randn(3)) ),
       "in_out"   => ( randn(srand(3),4), fft(randn(srand(3),3,4)*randn(srand(3),4) )  )
	     ),
	 Dict(
       "Operator" => ((*),),
       "params"   => (( 5*DFT((10,5)), -7*DCT((10,5))  ,),),
       "args"     => (  randn(10,5), fft(randn(10,5)) ),
       "in_out"   => ( randn(srand(3),10,5), 5*fft(-7*dct(randn(srand(3),10,5)))  )
	     ),
	 Dict(
       "Operator" => ((*),),
       "params"   => (( DFT((10,5)), -DCT((10,5)), GetIndex((10,10),(:,1:5))   ,),),
       "args"     => (  randn(10,10), fft(randn(10,5)) ),
       "in_out"   => ( randn(srand(3),10,10), fft(-dct(randn(srand(3),10,10)[:,1:5] ))  )
	     ),
	 Dict(
       "Operator" => ((*),),
       "params"   => (( DFT((10*5,)), Reshape(-DCT((10,10))[:,1:5],10*5)   ,),),
       "args"     => (  randn(10,10), fft(randn(10*5)) ),
       "in_out"   => ( randn(srand(3),10,10), fft(reshape(-dct(randn(srand(3),10,10))[:,1:5],10*5) )  )
	     ),
	 Dict(
       "Operator" => ((.*),),
       "params"   => (( randn(srand(1),3) , DCT((3,))                     ,),),
       "wrg_pr"   => (( randn(srand(1),5) , DCT((3,))                     ,),),
       "args"     => ( randn(3), randn(3) ),
       "in_out"   => ( randn(srand(3),3), randn(srand(1),3).*dct( randn(srand(3),3) )  )
	     ),
### testing HCAT ####
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( MatrixOp(randn(srand(1),3,10)), MatrixOp(randn(srand(2),3,4)),),),
       "wrg_pr"   => (( MatrixOp(randn(srand(1),5,10)), MatrixOp(randn(srand(2),3,4)),),),
       "args"     => ( [randn(10),randn(4)], randn(3) ),
       "in_out"   => ( [randn(srand(3),10),randn(srand(4),4)], 
		        randn(srand(1),3,10)*randn(srand(3),10)+
			randn(srand(2),3,4)*randn(srand(4),4)    )
	     ),
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( MatrixOp(randn(srand(1),3,10)), -MatrixOp(randn(srand(2),3,4)),),),
       "args"     => ( [randn(10),randn(4)], randn(3) ),
       "in_out"   => ( [randn(srand(3),10),randn(srand(4),4)], 
		        randn(srand(1),3,10)*randn(srand(3),10)-
			randn(srand(2),3,4)*randn(srand(4),4)    )
	     ),
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( DFT(6)', DCT(6),),),
       "wrg_pr"   => ((IDFT(6),  DCT(6),),),
       "args"     => ( Array[fft(randn(6)),randn(6)], randn(6) ),
       "in_out"   => ( Array[fft(randn(srand(1),6)),randn(srand(1),6)], 
		      real(6*ifft(fft(randn(srand(1),6)))+dct(randn(srand(1),6)))
		      ),
	     ),
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( MatrixOp(randn(srand(1),3,10)), -MatrixOp(randn(srand(2),3,4)), DCT(3) ),),
       "args"     => ( [randn(10),randn(4),randn(3)], randn(3) ),
       "in_out"   => ( [randn(srand(3),10),randn(srand(4),4),randn(srand(5),3)], 
		        randn(srand(1),3,10)*randn(srand(3),10)-
			randn(srand(2),3,4)*randn(srand(4),4)+
		       dct(randn(srand(5),3)))
	     ),
	 Dict(
       "Operator" => ((+),),
       "params"   => (( 
		       [MatrixOp(randn(srand(1),3,10)) -MatrixOp(randn(srand(2),3,4))],
		       [MatrixOp(randn(srand(3),3,10)) -MatrixOp(randn(srand(4),3,4))],
		       ),),
       "wrg_pr"   => (( 
		       [MatrixOp(randn(srand(1),3,10)) -MatrixOp(randn(srand(2),3,4))],
		       [MatrixOp(randn(srand(3),4,10)) -MatrixOp(randn(srand(4),4,4))],
		       ),),
       "args"     => ( [randn(10),randn(4)], randn(3) ),
       "in_out"   => ( [randn(srand(10),10),randn(srand(11),4)], 
		        randn(srand(1),3,10)*randn(srand(10),10)-
			randn(srand(2),3,4)*randn(srand(11),4)+
		        randn(srand(3),3,10)*randn(srand(10),10)-
			randn(srand(4),3,4)*randn(srand(11),4) )
	     ),
	 Dict(
       "Operator" => ((-),),
       "params"   => (( 
		       [MatrixOp(randn(srand(1),3,10)) -MatrixOp(randn(srand(2),3,4))],
		       [MatrixOp(randn(srand(3),3,10)) -MatrixOp(randn(srand(4),3,4))],
		       ),),
       "args"     => ( [randn(10),randn(4)], randn(3) ),
       "in_out"   => ( [randn(srand(10),10),randn(srand(11),4)], 
		        randn(srand(1),3,10)*randn(srand(10),10)-
			randn(srand(2),3,4)*randn(srand(11),4)-
		        randn(srand(3),3,10)*randn(srand(10),10)+
			randn(srand(4),3,4)*randn(srand(11),4) )
	     ),
## testing HCAT #### #TODO add some more tests on HCAT
	 Dict(
       "Operator" => ((vcat),),
       "params"   => (( MatrixOp(randn(srand(1),10,3)), MatrixOp(randn(srand(2),4,3)),),),
       "wrg_pr"   => (( MatrixOp(randn(srand(1),10,4)), MatrixOp(randn(srand(2),4,3)),),),
       "args"     => ( randn(3), [randn(10),randn(4)]),
       "in_out"   => ( randn(srand(3),3), 
		        [randn(srand(1),10,3)*randn(srand(3),3),
	                 randn(srand(2),4, 3)*randn(srand(3),3)]    )
	     ),
	 ]


for i in eachindex(stuff)

	x,y = deepcopy(stuff[i]["args"])

	params = stuff[i]["params"][1]
	Op     = stuff[i]["Operator"][1]
	A = Op(params...)

	test1,test2 = RegLS.test_FwAdj(A, x, y)
	@test test1 < 1e-8
	@test test2 < 1e-8
	test3 = RegLS.test_Op(A, x, y)
	@test test3 < 1e-8

	if "in_out" in keys(stuff[i]) 
		println("testing output")
		x0,y0 = stuff[i]["in_out"]
		test4 = norm(A*x0-y0)
		@test test4 < 1e-8
	end

	if "wrg_pr" in keys(stuff[i]) 
		params = stuff[i]["wrg_pr"][1]
		Op     = stuff[i]["Operator"][1]
		@test_throws Exception Op(params...)
	end

end

#
#@test norm(A*x-(A1.A*x-dct(x))) < 1e-8
#
#A1 = MatrixOp(randn(3,3))
#A2 = DCT(randn(3))
#A3 = MatrixOp(randn(3,3))
#
#A = A1-A2+A3
#
#test1,test2 = RegLS.test_FwAdj(A, x, y)
#@test test1 < 1e-8
#@test test2 < 1e-8
#test3 = RegLS.test_Op(A, x, y)
#@test test3 < 1e-8
#
#@test norm(A*x-(A1.A*x-dct(x)+A3.A*x   )) < 1e-8
#
#A = A1-A2-A3
#@test norm(A*x-(A1.A*x-dct(x)-A3.A*x   )) < 1e-8
#A = A1-(A2+A3)
#@test norm(A*x-(A1.A*x-dct(x)-A3.A*x   )) < 1e-8
#A = (A1-A2)-A3
#@test norm(A*x-(A1.A*x-dct(x)-A3.A*x   )) < 1e-8
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
