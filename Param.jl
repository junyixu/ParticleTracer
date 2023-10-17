module Param

# main(P::Parameters) 中会用到
export Parameters

# 必须要 using .Pushers 后
# Parameters 的 pusher 的类型才能是函数
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
