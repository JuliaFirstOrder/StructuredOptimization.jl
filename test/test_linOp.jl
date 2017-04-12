##test linear operators

stuff = [
## testing constructors ###
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
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => ((Float64,20,Float64,10, 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),20,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),20,10),b)
		      ),),
       "args"     => ( randn(10), randn(20) ),
	     ),
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => ((Float64,20,10, 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),20,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),20,10),b)
		      ),),
       "args"     => ( randn(10), randn(20) ),
	     ),
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => ((Float64,(20,),(10,), 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),20,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),20,10),b)
		      ),),
       "args"     => ( randn(10), randn(20) ),
	     ),
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => ((20,10, 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),20,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),20,10),b)
		      ),),
       "args"     => ( randn(10), randn(20) ),
	     ),
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => (((20,),(10,), 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),20,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),20,10),b)
		      ),),
       "args"     => ( randn(10), randn(20) ),
	     ),
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => ((10, 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),10,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),10,10),b)
		      ),),
       "args"     => ( randn(10), randn(10) ),
	     ),
	 Dict(
       "Operator" => (MyOperator,),
       "params"   => (((10,), 
		       (y,b) ->  A_mul_B!(y,randn(srand(1),10,10),b), 
		       (y,b) -> Ac_mul_B!(y,randn(srand(1),10,10),b)
		      ),),
       "args"     => ( randn(10), randn(10) ),
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

verb = true

for i in eachindex(stuff)

	x,y = deepcopy(stuff[i]["args"])

	params = stuff[i]["params"][1]
	Op     = stuff[i]["Operator"][1]
	A = Op(params...)

	test1,test2 = RegLS.test_FwAdj(A, x, y, verb)
	@test test1 < 1e-8
	@test test2 < 1e-8
	test3 = RegLS.test_Op(A, x, y)
	@test test3 < 1e-8

	if "in_out" in keys(stuff[i]) 
		verb && println("testing output")
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
