module PtcStruct

export Particle

struct Particle
	X::Matrix
	V::Matrix
	m::Float64
end

# function init_cond(p::Particle)
# 	echorow
# end

end
