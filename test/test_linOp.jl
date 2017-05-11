##test linear operators

stuff = [
#### testing constructors ###
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
       "Operator" => (Eye,),
       "params"   => ((randn(4,4),),),
       "args"     => ( randn(4,4), randn(4,4) ),
	     ),
	 Dict(
       "Operator" => (ZeroPad,),
       "params"   => (((4,),3),),
       "wrg_pr"   => (((4,),-3),),
       "args"     => ( randn(4), randn(7) ),
       "in_out"     => ( ones(4), [ones(4);zeros(3)] ),
	     ),
	 Dict(
       "Operator" => (ZeroPad,),
       "params"   => (((4,4),3),),
       "args"     => ( randn(4,4), randn(7,4) ),
       "in_out"     => ( ones(4,4), [ones(4,4);zeros(3,4)] ),
	     ),
	 Dict(
       "Operator" => (ZeroPad,),
       "params"   => (((4,4),3,3),),
       "args"     => ( randn(4,4), randn(7,7) ),
       "in_out"     => ( ones(4,4), [[ones(4,4);zeros(3,4)] zeros(7,3)] ),
	     ),
	 Dict(
       "Operator" => (ZeroPad,),
       "params"   => (((4,4,4),3,3),),
       "args"     => ( randn(4,4,4), randn(7,7,4) ),
       "in_out"     => ( ones(4,4,4), [[ones(4,4,4);zeros(3,4,4)] zeros(7,3,4)] ),
	     ),
	 Dict(
       "Operator" => (GetIndex,),
       "params"   => (( zeros(4,4), (1:4,)),),
       "args"     => ( randn(4,4), randn(4) ),
       "in_out"   => ( ones(4,4), ones(4) )
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
       "params"   => ((Complex{Float64}, (4,4), (1:2,)  ,),),
       "args"     => ( randn(4,4)+im, randn(2)+im ),
       "in_out"   => ( ones(4,4)+im, ones(2)+im )
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
       "Operator" => (Filt,),
       "params"   => ((300, [1.;0.;1.;0.;0.], [1.;1.;1.]),),
       "args"     => ( randn(300), randn(300) ),
       "in_out"   => ( randn(srand(1),300), filt([1.;0.;1.], [1.;1.;1.],randn(srand(1),300)) )
	     ),
	 Dict(
       "Operator" => (Filt,),
       "params"   => (((300,2), [1.;0.;1.;0.;0.], [1.;1.;1.]),),
       "args"     => ( randn(300,2), randn(300,2) ),
       "in_out"   => ( randn(srand(1),300,2), filt([1.;0.;1.], [1.;1.;1.],randn(srand(1),300,2)) )
	     ),
	 Dict(
       "Operator" => (Filt,),
       "params"   => ((300, [1.;0.;1.;0.;0.], [10.]),),
       "args"     => ( randn(300), randn(300) ),
       "in_out"   => ( randn(srand(1),300), filt([1.;0.;1.], [10.],randn(srand(1),300)) )
	     ),
	 Dict(
       "Operator" => (Filt,),
       "params"   => (((300,2), [1.;0.;1.;0.;0.], [10.]),),
       "args"     => ( randn(300,2), randn(300,2) ),
       "in_out"   => ( randn(srand(1),300,2), filt([1.;0.;1.], [10.],randn(srand(1),300,2)) )
	     ),
	 Dict(
       "Operator" => (MIMOFilt,),
       "params"   => (((10,3), [randn(10),randn(10),randn(10),randn(10),randn(10),randn(10)]),),
       "args"     => ( randn(10,3), randn(10,2) ),
	     ),

	 Dict(
       "Operator" => (MIMOFilt,),
       "params"   => (((10,2), [[1.;0.;1.;0.;0.],[1.;0.;1.;0.;0.]], 
		       [[1.;1.;1.],[2.;2.;2.]]),),
       "args"     => ( randn(10,2), randn(10,1) ),
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
## testing Reshape ####
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
### testing Scale ####
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
       "params"   => ((MatrixOp(2*ones(4,8)),),),
       "args"     => ( randn(8), randn(4) ),
       "in_out"   => ( ones(8), -2*ones(4,8)*ones(8) ),
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
       "Operator" => ((x -> broadcast(*, x... ),)),
       "params"   => (( (randn(srand(1),3) , FiniteDiff((3,)))                     ,),),
       "wrg_pr"   => (( (randn(srand(1),5) , FiniteDiff((3,)))                     ,),),
       "args"     => ( randn(3), randn(3) ),
       "in_out"   => ( randn(srand(3),3), randn(srand(1),3).*(FiniteDiff((3,))*randn(srand(3),3) )  )
	     ),
## testing Compose special cases ####
	 Dict(
       "Operator" => ((*),),
       "params"   => (( DFT((10*5,)), Reshape(-DCT((10,10))[:,1:5],10*5)   ,),),
       "args"     => (  randn(10,10), fft(randn(10*5)) ),
       "in_out"   => ( randn(srand(3),10,10), fft(reshape(-dct(randn(srand(3),10,10))[:,1:5],10*5) )  )
	     ),
### testing HCAT ####
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( MatrixOp(randn(srand(1),3,10)), MatrixOp(randn(srand(2),3,4)),),),
       "wrg_pr"   => (( MatrixOp(randn(srand(1),5,10)), MatrixOp(randn(srand(2),3,4)),),),
       "args"     => ( (randn(10),randn(4)), randn(3) ),
       "in_out"   => ( (randn(srand(3),10),randn(srand(4),4)), 
		        randn(srand(1),3,10)*randn(srand(3),10)+
			randn(srand(2),3,4)*randn(srand(4),4)    )
	     ),
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( MatrixOp(randn(srand(1),3,10)), -MatrixOp(randn(srand(2),3,4)),),),
       "args"     => ( (randn(10),randn(4)), randn(3) ),
       "in_out"   => ( (randn(srand(3),10),randn(srand(4),4)), 
		        randn(srand(1),3,10)*randn(srand(3),10)-
			randn(srand(2),3,4)*randn(srand(4),4)    )
	     ),
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( DFT(6)', DCT(6),),),
       "wrg_pr"   => ((IDFT(6),  DCT(6),),),
       "args"     => ( (fft(randn(6)),randn(6)), randn(6) ),
       "in_out"   => ( (fft(randn(srand(1),6)),randn(srand(1),6)), 
		      real(6*ifft(fft(randn(srand(1),6)))+dct(randn(srand(1),6)))
		      ),
	     ),
	 Dict(
       "Operator" => ((hcat),),
       "params"   => (( MatrixOp(randn(srand(1),3,10)), -MatrixOp(randn(srand(2),3,4)), DCT(3) ),),
       "args"     => ( (randn(10),randn(4),randn(3)), randn(3) ),
       "in_out"   => ( (randn(srand(3),10),randn(srand(4),4),randn(srand(5),3)), 
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
       "args"     => ( (randn(10),randn(4)), randn(3) ),
       "in_out"   => ( (randn(srand(10),10),randn(srand(11),4)), 
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
       "args"     => ( (randn(10),randn(4)), randn(3) ),
       "in_out"   => ( (randn(srand(10),10),randn(srand(11),4)), 
		        randn(srand(1),3,10)*randn(srand(10),10)-
			randn(srand(2),3,4)*randn(srand(11),4)-
		        randn(srand(3),3,10)*randn(srand(10),10)+
			randn(srand(4),3,4)*randn(srand(11),4) )
	     ),
## testing VCAT #### #TODO add some more tests on VCAT
	 Dict(
       "Operator" => ((vcat),),
       "params"   => (( MatrixOp(randn(srand(1),10,3)), MatrixOp(randn(srand(2),4,3)),),),
       "wrg_pr"   => (( MatrixOp(randn(srand(1),10,4)), MatrixOp(randn(srand(2),4,3)),),),
       "args"     => ( randn(3), (randn(10),randn(4))),
       "in_out"   => ( randn(srand(3),3), 
		        (randn(srand(1),10,3)*randn(srand(3),3),
	   randn(srand(2),4, 3)*randn(srand(3),3))    )
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
		test4 = RegLS.deepvecnorm(A*x0.-y0)
		@test test4 < 1e-8
	end

	if "wrg_pr" in keys(stuff[i]) 
		params = stuff[i]["wrg_pr"][1]
		Op     = stuff[i]["Operator"][1]
		@test_throws Exception Op(params...)
	end

end

##testing some special cases...

L = 3*DiagOp(3*ones(3))
@test typeof(L) <: DiagOp{Vector{Float64},1}
@test norm(L.d-9)<1e-9

L = (3+3*im)*DiagOp(3*ones(3))
@test typeof(L) <: DiagOp{Vector{Complex{Float64}},1}
@test norm(L.d-(9+im*9))<1e-9

L = (3*ones(3)).*(3*Eye(3))
@test typeof(L) <: DiagOp{Vector{Float64},1}
@test norm(L.d-(9))<1e-9

L = DiagOp(3*ones(3))*(3*Eye(3))
@test typeof(L) <: DiagOp{Vector{Float64},1}
@test norm(L.d-(9))<1e-9



#TODO move this up

B = [randn(10),randn(5),randn(10),randn(2),randn(10),randn(10)]
A = [[1.],[1.],[1.],[1.],[1.],[1.]]

x = randn(100,3)
L = MIMOFilt((100,3), B )
y = L*x

@test norm(y[:,1]-(filt(B[1],A[1],x[:,1])+filt(B[2],A[2],x[:,2])+filt(B[3],A[3],x[:,3])   )) <1e-8
@test norm(y[:,2]-(filt(B[4],A[4],x[:,1])+filt(B[5],A[5],x[:,2])+filt(B[6],A[6],x[:,3])   )) <1e-8


#A = [FiniteDiff((4,),1) MatrixOp(randn(4,4))]
#
#println(A)
#
#y = zeros(4)
#x1,x2 = randn(4),randn(4)
#
#A_mul_B!(y,A,(x1,x2))
#
#At = A'
#println(typeof(At))
#
#@time A_mul_B!(y,A,(x1,x2))
#@time A_mul_B!(y,A,(x1,x2))
#@time A_mul_B!(y,A,(x1,x2))
#@code_warntype A_mul_B!(y,A,(x1,x2))







