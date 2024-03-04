module Fields
include("Constants.jl")
using .Constants

function B_rϕz2xyz(BR::Float64,Bϕ::Float64,Bz::Float64, cosϕ::Float64, sinϕ::Float64)
    B1 = BR*cosϕ-Bϕ*sinϕ;
    B2 = BR*sinϕ+Bϕ*cosϕ;
    B3 = Bz;
    return [B1, B2, B3]
end


"""
环坐标下 tokamak 场
"""
function tokamak(x::Float64,y::Float64,z::Float64,q)
    # q: 安全因子

    R0 = Constants.R0

	R²=x^2 + y^2
	R=sqrt(R²)
    B1=(-q*R0*y+x*z)/(R²*q)
    B2=(q*R0*x+y*z)/(R²*q)
    B3=(-1+R0/R)/q

    return [B1, B2, B3]
end
end
