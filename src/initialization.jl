using Distributions
function func_para(x::T, m::T) where T<:AbstractFloat
    return 4x*(m^2-x^2)/m^4
end

function TORD2CYLD(r::T, θ::T, ζ::T,R0::T) where T<:AbstractFloat #环坐标转为柱坐标
    R = R0 + r * cos(θ)
    ϕ = -ζ
    z = r * sin(θ)
    return [R, ϕ, z]
end

function TORD2CYLD(rθζ::AbstractVector{T},R0::T) where T<:AbstractFloat #环坐标转为柱坐标
    TORD2CYLD(rθζ[1], rθζ[2], rθζ[3], R0)
end


function CYLD2CART(R::T, ϕ::T, z::T) where T<:AbstractFloat           # 柱坐标(r,  Φ,   z)转为笛卡尔坐标(x,  y, z)
	x= R * cos(ϕ)
	y = R * sin(ϕ)
    return [x, y, z]
end

function CYLD2CART(Rϕz::AbstractVector{<:AbstractFloat})          # 柱坐标(r,  Φ,   z)转为笛卡尔坐标(x,  y, z)
    CYLD2CART(Rϕz[1], Rϕz[2], Rϕz[3])
end

function TORD2CART(rθζ::AbstractVector{T}, R0::T) where T<:AbstractFloat # 环坐标转为笛卡尔坐标
	Rϕz=TORD2CYLD(rθζ, R0)
	CYLD2CART(Rϕz)
end

function GenRandNum_Para(m::AbstractFloat)
    while true
        x=m*rand()
        y=func_para(x,m)
        dScope=func_para(sqrt(1/3)*m,m)*rand()
    	dScope>y || return x
    end
end


function local_frame(B::AbstractVector{<:AbstractFloat})
    local_frame(Bx,By,Bz)
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
end

function SetParticleMomentum_Gyrocenter(x::T, y::T, z::T, px::T, py::T, pz::T, B::AbstractVector{<:AbstractFloat}) where T <: AbstractFloat
	γ = rand(Normal(μ=UserInputs.μ, σ=UserInputs.σ))                         # 产生一个服从正态分布的 gamma，范围在(Mu-5*Sigma,Mu+5*Sigma)
    p = sqrt(γ^2-1) # ||p||
	θ_min	= UserInputs.θ_min # 投掷角p⊥/p的最小值
	θ_max= UserInputs.θ_max # 投掷角p⊥/p的最大值
	ϕ_min	= UserInputs.ϕ_min # 角度的最小值
	ϕ_max	= UserInputs.ϕ_max # 角度的最大值
    ϕ = ϕ_min+(ϕ_max-ϕ_min)*rand()
    θ = θ_min+(θ_max-θ_min)*rand()
    sinθ = sin(θ)
    cosθ = sqrt(1 - sinθ^2)
	p⊥ = p*sinθ
    p_para = p*cosθ
	e1,e2,b = local_frame(B)
    return @. p_para * b + p⊥ * (cos(ϕ) * e1 + sin(ϕ) * e2)
    # Transform (p_\parallel,p_\perp,phase) to (px,py,pz)
end

