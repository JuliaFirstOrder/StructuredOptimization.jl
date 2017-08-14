export problem

function problem(terms::Vararg)
	cf = ()
	for i = 1:length(terms)
		cf = (cf...,terms[i]...)
	end
	return cf
end

# problem = vcat # why tuples? with arrays this is so much simpler
