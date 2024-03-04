module Pushers
using Test
using LinearAlgebra:norm
export pusher2, pusher1
using LinearAlgebra: I, ⋅ 

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

q=1
m=1
c=1
dt=0.1

function _Omega(B::AbstractVector)
	mat = [  0    -B[3]   B[2]
		    B[3]   0     -B[1]
		   -B[2]   B[1]    0   ]
	return q*dt/(2m*c) .* mat
end

function boris(x::AbstractVector, p::AbstractVector, E::AbstractVector, B::AbstractVector)
	Ω=_Omega(B)
	R=inv(I + Ω)*(I - Ω)
    γ = sqrt(1 + p⋅p/(m^2*c^2))
    # println(γ)
    # γ=1
    v = p/(γ*m)
	new_v = R*v + q*dt/m * inv(I+Ω) * E
	new_x = x + new_v*dt
	return (new_x,new_v*γ*m)
end
end
