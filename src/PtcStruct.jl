module PtcStruct
include("Constants.jl")
export Particle, ParticleData

struct Particle{M<:AbstractVector{<:AbstractFloat}}
	X::M
	P::M
	B::M
end

struct ParticleData{M<:AbstractMatrix{<:AbstractFloat}}
	X::M
	P::M
	B::M
end

end
