using Distributions
R0 = UserInputs.R0

function local_frame(B::AbstractVector{<:AbstractFloat})
    local_frame(B[1],B[2],B[3]) # Bx,By,Bz
end

function local_frame(Bx::T,By::T,Bz::T) where T<:AbstractFloat
    normB = sqrt(Bx^2 + By^2 + Bz^2) # 计算 磁场 B 的大小
    normB_xy= sqrt(Bx^2 + By^2) # 计算 磁场 B 在xoy面投影  的大小
	if normB==0
		error("In function local_frame: There is no magnetic field.")
    elseif normB_xy==0
        e1 = [1, 0, 0]
        e2 = [0, 1, 0]
        b = [0, 0, 1]
	else
        e1 = [Bz*Bx/(normB*normB_xy) 
              Bz*By/(normB*normB_xy) 
              -1*normB_xy/normB]
        e2 = [-1*By/normB_xy
              Bx/normB_xy
              0]
        b  = [Bx/normB
              By/normB
              Bz/normB]
    end
    return (e1,e2,b)
end

function SetParticleMomentum_Gyrocenter(x::T, y::T, z::T, B::AbstractVector{<:AbstractFloat}) where T <: AbstractFloat
	γ = rand(Distributions.Normal(UserInputs.μ, UserInputs.σ))                         # 产生一个服从正态分布的 gamma，范围在(Mu-5*Sigma,Mu+5*Sigma)
    if γ^2-1 < 0 
        (γ = UserInputs.μ)
        @warn("γ^2-1 < 0!")
    end
    p = sqrt(γ^2-1) # ||p||
	sinθ_min	= UserInputs.sinθ_min # p⊥/p的最小值
	sinθ_max= UserInputs.sinθ_max # p⊥/p的最大值
	ϕ_min	= UserInputs.ϕ_min
	ϕ_max	= UserInputs.ϕ_max
    ϕ = ϕ_min+(ϕ_max-ϕ_min)*rand()
    sinθ = sinθ_min+(sinθ_max-sinθ_min)*rand()
    cosθ = sqrt(1 - sinθ^2)
	p⊥ = p*sinθ
    p_para = p*cosθ
	e1,e2,b = local_frame(B)
    return @. p_para * b + p⊥ * (cos(ϕ) * e1 + sin(ϕ) * e2)
    # Transform (p_\parallel,p_\perp,phase) to (px,py,pz)
end

