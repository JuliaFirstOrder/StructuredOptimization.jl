type LeastSquares{T <: AffineOp} <: SmoothTerm
	A::T
end

ls{T <: AffineOp}(A::T) = LeastSquares(A)
ls(x::OptVar) = LeastSquares(eye(x))

export ls
