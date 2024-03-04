module PtcStruct
include("Constants.jl")
export Particle

struct Particle{M<:AbstractMatrix}
	X::M
	P::M
	B::M
end

end
