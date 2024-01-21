module Pushers
using Test
using LinearAlgebra:norm
export pusher2, pusher1

include("parameters.jl")
const G = 1

function pusher_gravity(ps::Vector, i::Int, s::Int)
	# i for i-th particle
	# s for s-th time step
	xᵢ = ps[i].X[s,:]
	@test xᵢ isa Vector
	Σ = zeros(2)
	for j = [j for j = 1:length(ps) if j != i] # other particles except i
		mⱼ = ps[j].m
		xⱼ = ps[j].X[s,:]
		Σ += mⱼ * (xⱼ - xᵢ) / norm(xⱼ - xᵢ)^3
	end
	a = G * Σ

	vᵢ = ps[i].V[s,:]

	v_ = a * Δt + vᵢ
	x_ = a * Δt^2 + vᵢ * Δt + xᵢ
	return (x_, v_)
end

function pusher1(x::T, p::T) where T <:AbstractFloat
	x_ = x + p*Δt
	p_ = p - 0.5ω2*(x_ + x) * Δt
	return (x_,p_)
end

function pusher2(x::T, p::T) where T <:AbstractFloat
	x_ = x + p*Δt
	p_ = p - ω2*x_ * Δt
	return (x_,p_)
end

end
