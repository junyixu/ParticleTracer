module Pushers
export pusher2, pusher1

include("parameters.jl")

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
