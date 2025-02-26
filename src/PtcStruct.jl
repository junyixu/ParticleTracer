module PtcStruct
include("Constants.jl")
include("UserInputs.jl")
export Particle, ParticleData, MagneticParticle, EMParticle, AbstractParticle, AbstractParticleData, MagneticParticleData, EMParticleData

"""
抽象粒子类型，作为所有具体粒子类型的父类
"""
abstract type AbstractParticle end

"""
只考虑磁场的带电粒子
"""
struct MagneticParticle <: AbstractParticle
	X::Vector{Float64}  # 位置
	P::Vector{Float64}  # 动量
	B::Vector{Float64}  # 磁场
end

"""
同时考虑电磁场的带电粒子
"""
struct EMParticle <: AbstractParticle
	X::Vector{Float64}  # 位置
	P::Vector{Float64}  # 动量
	B::Vector{Float64}  # 磁场
	E::Vector{Float64}  # 电场
end

# 构造函数
function Particle(x0::AbstractVector, p0::AbstractVector, B0::AbstractVector)
	MagneticParticle(Vector{Float64}(x0), Vector{Float64}(p0), Vector{Float64}(B0))
end

function Particle(x0::AbstractVector, p0::AbstractVector, B0::AbstractVector, E0::AbstractVector)
	EMParticle(Vector{Float64}(x0), Vector{Float64}(p0), Vector{Float64}(B0), Vector{Float64}(E0))
end

"""
抽象粒子数据类型
"""
abstract type AbstractParticleData end

"""
只考虑磁场的粒子数据
"""
struct MagneticParticleData <: AbstractParticleData
	X::Matrix{Float64}
	P::Matrix{Float64}
	B::Matrix{Float64}
end

"""
同时考虑电磁场的粒子数据
"""
struct EMParticleData <: AbstractParticleData
	X::Matrix{Float64}
	P::Matrix{Float64}
	B::Matrix{Float64}
	E::Matrix{Float64}
end

# ParticleData 的构造函数
function ParticleData(x::AbstractMatrix, p::AbstractMatrix, b::AbstractMatrix)
	MagneticParticleData(Matrix{Float64}(x), Matrix{Float64}(p), Matrix{Float64}(b))
end

function ParticleData(x::AbstractMatrix, p::AbstractMatrix, b::AbstractMatrix, e::AbstractMatrix)
	EMParticleData(Matrix{Float64}(x), Matrix{Float64}(p), Matrix{Float64}(b), Matrix{Float64}(e))
end

end
