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

"""
生成服从参数为m的抛物线分布的随机数
m 一般为小半径 r_max = a
"""
function Para(m::AbstractFloat) 
    while true
        x=m*rand()
        y=func_para(x,m)
        dScope=func_para(sqrt(1/3)*m,m)*rand()
    	dScope>y || return x
    end
end


function SetParticlePosition_ParabolicTorus(r_max::AbstractFloat)
     rθζ=[Para(r_max)
          2π*rand()
          2π*rand()]
	TORD2CART(rθζ,R0)
end
