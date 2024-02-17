module PtcStruct

export Particle

struct Particle{M<:AbstractMatrix}
	X::M
	V::M
	m::Float64
end

# function init_cond(p::Particle)
# 	echorow
# end

end
