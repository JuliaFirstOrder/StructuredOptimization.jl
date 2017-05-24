@printf("\nTesting linear operators\n")

function test_op(A::LinearOperator, x, y, verb::Bool = false)
	verb && (println(); show(A); println())

  Ax = A*x
	Ax2 = RegLS.deepsimilar(Ax)
  verb && println("forward preallocated")
  A_mul_B!(Ax2, A, x) #verify in-place linear operator works
  verb && @time A_mul_B!(Ax2, A, x)

  @test RegLS.deepvecnorm(Ax .- Ax2) <= 1e-8

  Acy = A'*y
	Acy2 = RegLS.deepsimilar(Acy)
  verb && println("adjoint preallocated")
  Ac_mul_B!(Acy2, A, y) #verify in-place linear operator works
  verb && @time Ac_mul_B!(Acy2, A, y)

  @test RegLS.deepvecnorm(Acy .- Acy2) <= 1e-8

	s1 = RegLS.deepvecdot(Ax2, y)
	s2 = RegLS.deepvecdot(x, Acy2)

	@test abs( s1 - s2 ) < 1e-8

  return Ax
end

# Conv

# DCT

# DFT

# DiagOp

# Eye

# Filt

# FiniteDiff

# GetIndex

# MatrixOp

# MIMOFilt

# Variation

# Xcorr

# ZeroPad

# Zeros
