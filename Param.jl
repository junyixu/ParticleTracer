module Param

export Parameters

include("Pushers.jl")
using .Pushers
mutable struct Parameters{T<:AbstractFloat}
	filename::String
	nsteps::Int
	nfiles::Int
	ω2::T
	x0::T
	p0::T
	Δt::T
	pusher::Function
end

include("parameters.jl")

p=Parameters(filename,
			 5000,
			 20,
			 ω2,
			 x0,
			 p0,
			 Δt,
			 @eval $pusher)
end
