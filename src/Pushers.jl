module Pushers
using Test
using LinearAlgebra:norm
export pusher2, pusher1
using LinearAlgebra: I, ⋅ 

include("UserInputs.jl")
using .UserInputs: Δt

q=1
m=1
c=1

function _Omega(B::AbstractVector)
	mat = [  0    -B[3]   B[2]
		    B[3]   0     -B[1]
		   -B[2]   B[1]    0   ]
	return q*Δt/(2m*c) .* mat
end

function boris(x::AbstractVector, p::AbstractVector, E::AbstractVector, B::AbstractVector)
	Ω=_Omega(B)
	R=inv(I + Ω)*(I - Ω)
    γ = sqrt(1 + p⋅p/(m^2*c^2))
    # println(γ)
    # γ=1
    v = p/(γ*m)
	new_v = R*v + q*Δt/m * inv(I+Ω) * E
	new_x = x + new_v*Δt
	return (new_x,new_v*γ*m)
end

function p₋2p₊(B::Vector{T}, γ::T, p::Vector{T})::Vector{T} where T <:AbstractFloat
    # B[1]: Bx
	# B[2]: By
	# B[3]: Bz
	# p[1]: px
	# p[2]: py
	# p[3]: pz
    Bx = B[1]
	By = B[2]
	Bz = B[3]
	px = p[1]
	py = p[2]
	pz = p[3]

	temp  = Bx^2*Δt^2 + By^2*Δt^2 + Bz^2*Δt^2 + 4*γ^2

	px₊ = (Bx^2*Δt^2*px - By^2*Δt^2*px - Bz^2*Δt^2*px + 2*Bx*By*Δt^2*py +  2*Bx*Bz*Δt^2*pz + 4*Bz*Δt*py*γ - 4*By*Δt*pz*γ + 4*px*γ^2)/temp

	py₊ =(2*Bx*By*Δt^2*px - Bx^2*Δt^2*py + By^2*Δt^2*py - Bz^2*Δt^2*py +  2*By*Bz*Δt^2*pz - 4*Bz*Δt*px*γ + 4*Bx*Δt*pz*γ + 4*py*γ^2)/temp

	pz₊ =(2*Bx*Bz*Δt^2*px + 2*By*Bz*Δt^2*py - Bx^2*Δt^2*pz - By^2*Δt^2*pz +   Bz^2*Δt^2*pz + 4*By*Δt*px*γ - 4*Bx*Δt*py*γ + 4*pz*γ^2)/temp

    return [px₊, py₊, pz₊] # new p
end

function RVPA_Cay3D(x::AbstractVector{T}, p::AbstractVector{T}, E::AbstractVector{T}, B::AbstractVector{T}) where T <: AbstractFloat
	# p₋
    p₋ = p + 0.5Δt*E;

    γ = 1 + p₋⋅p₋ # m = 1; c = 1

    p₊ = p₋2p₊(B, γ, p₋)

    p₊ = p₊ + 0.5Δt*E;

    γ = 1 + p₊⋅p₊ # m = 1; c = 1

    v₊ =  p₊ / γ
    x₊ = x + v₊ * Δt

    return (x₊, p₊)
end

end
