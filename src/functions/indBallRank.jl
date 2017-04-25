immutable IndBallRank{I <: Integer}   <: NonSmoothFunction
	r::I
end

get_prox(T::IndBallRank) = ProximalOperators.IndBallRank(T.r)

fun_name(T::IndBallRank, i::Int64) = "Ind{rank(A$(i)x) ≤ r$(i)}(x) "
fun_par( T::IndBallRank, i::Int64) = " r$(i) = $(T.r) "
