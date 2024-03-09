module PtcStruct
include("Constants.jl")
export Particle

struct Particle{M<:AbstractMatrix{<:AbstractFloat}}
	X::M
	P::M
	B::M
end

end
