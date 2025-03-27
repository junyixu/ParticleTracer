#=
    Data analysis for APT
    Copyright © 2023 Junyi Xu <jyxu@mail.ustc.edu.cn>

    Distributed under terms of the MIT license.
=#

module Ptcs
using LinearAlgebra:⋅, × 
using HDF5
using OffsetArrays
include("InTriangular.jl")
using .InTriangular
const mα=6.6446573357e-27
const me=9.1093837015e-31

const τA_E = 6.59e-7
const τA_I = 6.35e-7

import Base: size, length
export Unit, Particle, Particles, MetaData, is_in_flux_surface, get_psi_by_pos, get_vertex_id_by_pos, cart2cyld, cart2polar, is_dead, is_passage, ParticleEB, ParticleOnlyB

struct Unit
    B::Float64
    E::Float64 # 电场
    Ω::Float64
    t::Float64
    x::Float64 # 位置
    p::Float64
    v::Float64
    ε::Float64 # 能量
    q::Float64
    m::Float64
    function Unit(B::Float64, q::Float64, m::Float64)
		c_0 = 2.99792458e8
        E = B*c_0;
        Ω=q*B/m;
        t= 1/Ω;
        x = m*c_0/(q*B);
        p=m*c_0;
        v=c_0;
        ε=m*c_0^2;
        new(B, E, Ω, t, x, p, v, ε, q, m)
    end
end

function Unit(B::Float64, type::String)
	if type == "alpha"
		return Unit(B, 2*1.60217733e-19, mα) # 有时候需要改单位电量和单位质量
	end
	return Unit(B, 1.60217733e-19, me)
end
Unit(B::Float64)=Unit(B, 1.60217733e-19, 9.1093837015e-31)

abstract type Particle end

struct ParticleEB <: Particle
	dT::Float64
	die_step::Int
	isDead::Bool
	X::Matrix
	P::Matrix
	B::Matrix
	E::Matrix
	l::Int # 总步数
	function ParticleEB(PTC::Matrix, B::Matrix, E::Matrix, u::Unit)
		die = PTC[:, 1] .!= 0.0
		dT = PTC[2,2] * u.t
		die_step = 1
		isDead = false
		for i in eachindex(die) # 在第几步死亡
			if die[i] == 1.0
				die_step = i
				isDead = true
				break
			end
		end
		if !isDead
			die_step = length(die)+1
		end
		X = PTC[:, 3:5] * u.x
		P = PTC[:, 6:8] * u.p
		new(dT, die_step, isDead, X, P, B*u.B, E*u.E, size(PTC, 1))
	end
end

struct ParticleOnlyB <: Particle
	dT::Float64
	die_step::Int
	isDead::Bool
	X::Matrix
	P::Matrix
	B::Matrix
	l::Int # 总步数
	function ParticleOnlyB(PTC::Matrix, B::Matrix, u::Unit)
		die = PTC[:, 1] .!= 0.0
		dT = PTC[2,2] * u.t
		die_step = 1
		isDead = false
		for i in eachindex(die) # 在第几步死亡
			if die[i] == 1.0
				die_step = i
				isDead = true
				break
			end
		end
		if !isDead
			die_step = length(die)+1
		end
		X = PTC[:, 3:5] * u.x
		P = PTC[:, 6:8] * u.p
		new(dT, die_step, isDead, X, P, B*u.B, size(PTC, 1))
	end
end

function Particle(PTC::Matrix, B::Matrix, u::Unit)
    ParticleOnlyB(PTC, B, u)
end

function Particle(PTC::Matrix, B::Matrix, E::Matrix, u::Unit)
    ParticleEB(PTC, B, E, u)
end

function Particle(filename::String, u::Unit)
	fid=h5open(filename, "r")
	PTC = read(fid["PTC"])
	B = read(fid["B"])
    if haskey(fid, "E")
        E = read(fid["E"])
        close(fid)
        return ParticleEB(PTC, B, E, u)
    else
        close(fid)
        return ParticleOnlyB(PTC, B, u)
    end
end

# 初始 pitch angle
function cos(P::Matrix, B::Matrix)::Vector # P: 3xn B: 3xn
	# 初始磁场
	B_norm = sqrt.(sum(B.^2, dims=1))
	# 初始动量
	p_norm = sqrt.(sum(P.^2, dims=1))
	# 初始 pitch angle
	return sum(P.*B, dims=1) ./ B_norm ./ p_norm |> vec
end

function cos(p::Vector, B::Vector)::AbstractFloat
	B_norm = (sqrt ∘ sum)(B.^2)
	p_norm = (sqrt ∘ sum)(p.^2)
	return sum(p.*B) / B_norm / p_norm
end

# %%
struct Particles
	ptcs::Vector{Particle}
	dT::Float64 # 步长/秒 
	X0::Matrix # 初始位置/米
	P0::Matrix # 初始动量
	B0::Matrix # 初始磁场/特斯拉
	die::Vector # 在第i步死亡的粒子数
	cosθ::Vector 
	n::Int # 粒子数
	l::Int # 总步数
	function Particles(filename::String, u::Unit)
		fid=h5open(filename, "r")
		# for key in keys(fid) # 读取全部数据
		# 	sb=Symbol(key)
		# 	eval(:($sb = read(fid[$key])))
		# end
		PTC = read(fid["PTC"])
		B = read(fid["B"])
        if haskey(fid, "E")
            E = read(fid["E"])
        end

		dT = PTC[2,2] * u.t

		n = Int(size(PTC, 2)/11) # size, n ptcs
		l = size(PTC, 1) # size, n ptcs

		ptcs = Vector{Particle}(undef, n)
		for i in 1:Int(size(PTC, 2)/11)
            if haskey(fid, "E")
                ptcs[i] = Particle(PTC[:, (i-1)*11+1:i*11], B[:, (i-1)*3+1:i*3], E[:, (i-1)*3+1:i*3], u)
            else
                ptcs[i] = Particle(PTC[:, (i-1)*11+1:i*11], B[:, (i-1)*3+1:i*3], u)
            end
		end

		close(fid)
		die = zeros(Int, l)
		for i in 1:n
			if ptcs[i].isDead
				die[ptcs[i].die_step]+=1
			end
		end

		X0 = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.X[1, :])
		P0 = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.P[1, :])
		B0 = reshape(B[1, :], 3, :) * u.B
		new(ptcs, dT, X0, P0, B0, die, cos(P0, B0), n, l)
	end
	function Particles(ptcs::Vector{Particle})
		n = length(ptcs)
		l = ptcs[1].l
		die = zeros(Int, l)
		for i in 1:n
			if ptcs[i].isDead
				die[ptcs[i].die_step]+=1
			end
		end
		X0 = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.X[1, :])
		P0 = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.P[1, :])
		B0 = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.B[1, :])
		dT = ptcs[1].dT
		new(ptcs, dT, X0, P0, B0, die, cos(P0, B0), n, l)
	end
end
# %%

Base.length(ptc::Particle) = size(ptc.X , 1)
Base.size(ptc::Particle) = (size(ptc.X , 1),)

import Base: setindex!, getindex, broadcast, iterate, lastindex, length
setindex!(p::Particles, i) = setindex!(p.ptcs, i)
getindex(p::Particles, i) = getindex(p.ptcs, i)
broadcast(f::Function, p::Particles) = broadcast(f, p.ptcs)
lastindex(p::Particles) = lastindex(p.ptcs)
length(p::Particles) = length(p.ptcs)
iterate(p::Particles) = iterate(p.ptcs)
iterate(p::Particles, i) = iterate(p.ptcs, i)


x(ptc::Particle)::Vector = ptc.X[:, 1]
y(ptc::Particle)::Vector = ptc.X[:, 2]
z(ptc::Particle)::Vector = ptc.X[:, 3]
R(ptc::Particle)::Vector = sqrt.(x(ptc).^2 + y(ptc).^2)
function R(X::Matrix)
	if size(X, 2) != 3
		error("不是 (x,y,z) 矩阵")
	end
	@. sqrt(X[:, 1]^2 + X[:, 2])
end
function z(X::Matrix)
	if size(X, 2) != 3
		error("不是 (x,y,z) 矩阵")
	end
	X[:, 3]
end
px(ptc::Particle)::Vector = ptc.P[:, 1]
py(ptc::Particle)::Vector = ptc.P[:, 2]
pz(ptc::Particle)::Vector = ptc.P[:, 3]
p_norm(ptc::Particle)::Vector = sqrt.(px(ptc).^2 + py(ptc).^2 + pz(ptc).^2)

x(ptcs::Particles)::Matrix = hcat([x(ptcs.ptcs[i]) for i in 1:size(ptcs)]...)
y(ptcs::Particles)::Matrix = hcat([y(ptcs.ptcs[i]) for i in 1:size(ptcs)]...)
z(ptcs::Particles)::Matrix = hcat([z(ptcs.ptcs[i]) for i in 1:size(ptcs)]...)
R(ptcs::Particles)::Matrix = hcat([R(ptcs.ptcs[i]) for i in 1:size(ptcs)]...)
z(ptcs::Vector{Particle})::Matrix = hcat([z(ptcs[i]) for i in 1:length(ptcs)]...)
R(ptcs::Vector{Particle})::Matrix = hcat([R(ptcs[i]) for i in 1:length(ptcs)]...)

X(ptcs::Particles, i::Int) = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.X[i, :])
P(ptcs::Particles, i::Int) = reduce(hcat, [ptcs[i] for i in 1:n] .|> ptc->ptc.P[i, :])

is_dead(p::Particle)::Bool = p.isDead
is_dead(p::Particles)::BitVector = [is_dead(ptc) for ptc in p.ptcs]

struct MetaData
	T0::OffsetArray{Int}
	Tensor::OffsetArray{Int}
	T1::OffsetArray{Int}
	T2::OffsetArray{Int}
	XY::OffsetArray{Float64}
	d::Vector{Float64}
	o::Vector{Float64}
	s::Vector{Float64}
	function MetaData(filename::String)
		h5open(filename, "r") do fid
			T0= read(fid["T0/data"])
			Tensor = similar(T0, reverse(size(T0))) # 优化内存布局
			for i = 1:size(T0, 3)
				Tensor[i, :, :] = T0[:, :, i]'
			end
			Tensor=OffsetArray(Tensor, -1, -1, -1)
			T0=OffsetArray(T0, -1, -1, -1)
			T1=OffsetArray(read(fid["T1/data"]), -1, -1)
			T2=OffsetArray(read(fid["T2/data"]), -1, -1)
			XY=OffsetArray(read(fid["XY/data"]), -1, -1)
			new(T0, Tensor, T1, T2, XY,read(fid["d"]), read(fid["o"]), read(fid["s"]))
		end
	end
end


"""
通过 vertex 编号知道该 vertex 所对应的磁面 ψ 是第几层，
其中 20201 是 XY 总的点的个数，4 是公差
"""
i2ψ(i::Int)::Int=(i-1)%Int((20201 -1)/4) |> y -> 0.5*(-1 + sqrt(1+8y)) |> x->floor(Int, x)
# i2ψ(i::Int)::Int=(i-1)%Int((78805 -1)/4) |> y -> 0.5*(-1 + sqrt(1+8y)) |> x->floor(Int, x)

"""
判断点的坐标是否在最外层磁面内

	is_in_flux_surface(xy::AbstractVector, md::MetaData)

"""
function is_in_flux_surface(xy::AbstractVector, MD::MetaData) # 出 ITER 非截面
	i,j= floor.(Int, (parent(xy) .- MD.o)./MD.d)
	if 0<i<MD.s[1] && 0<j<MD.s[2]
		for vertex in MD.T0[i, j, 1:MD.T0[i, j, 0]]
			for tri in MD.T1[0:MD.T1[6, vertex]-1, vertex]
				if is_inside(parent(xy), parent(MD.XY[:, MD.T2[:, tri]]))
					return true
				end
			end
		end
	end
	return false
end

"""
通过点的坐标找到周围的数据点的编号

# Example

```julia-repl
julia> get_vertex_id_by_pos([r-R0,z]./u.x, MD)
3-element Vector{Int64}:
 14720
 14721
 14817
```
"""
function get_vertex_id_by_pos(xy::AbstractVector, MD::MetaData)::Vector{Int}
	i,j= floor.(Int, (parent(xy) .- MD.o)./MD.d)
	if 0<i<MD.s[1] && 0<j<MD.s[2]
		for vertex in MD.T0[i, j, 1:MD.T0[i, j, 0]]
            println(vertex)
			for tri in MD.T1[0:MD.T1[6, vertex]-1, vertex]
				if is_inside(parent(xy), parent(MD.XY[:, MD.T2[:, tri]]))
					return parent(MD.T2[:, tri])
				end
			end
		end
	end
	# error("出界!!!")
	return [-1,-1,-1]
end

function get_vertex_id_by_pos!(three::Vector{Int}, xy::Tuple{Float64, Float64}, MD::MetaData)::Int
	i,j= floor.(Int, (xy .- MD.o)./MD.d)
	if 0<i<MD.s[1] && 0<j<MD.s[2]
		for vertex in MD.T0[i, j, 1:MD.T0[i, j, 0]]
            if parent(MD.XY[:, vertex]) == xy
                return vertex
            end
			for tri in MD.T1[0:MD.T1[6, vertex]-1, vertex]
				if is_inside(xy, parent(MD.XY[:, MD.T2[:, tri]]))
					three .= parent(MD.T2[:, tri])
                    return -1
				end
			end
		end
	end
	# error("出界!!!")
	return -2
end

function get_vertex_id_by_pos!(three::Vector{Int}, xy::AbstractVector, MD::MetaData)::Int
	i,j= floor.(Int, (parent(xy) .- MD.o)./MD.d)
	# println("i = $i, j = $j")
	if 0<i<MD.s[1] && 0<j<MD.s[2]
		for vertex in MD.T0[i, j, 1:MD.T0[i, j, 0]]
            if parent(MD.XY[:, vertex]) == xy
                return vertex
            end
			for tri in MD.T1[0:MD.T1[6, vertex]-1, vertex]
				if is_inside(parent(xy), parent(MD.XY[:, MD.T2[:, tri]]))
					three .= parent(MD.T2[:, tri])
                    return -1
				end
			end
		end
        # return MD.T0[i, j, 1]
	end
	# error("出界!!!")
	return -2
end

function get_vertex_id_by_pos_type5(xy::AbstractVector, MD::MetaData)::Vector{Int}
    i,j = floor.(Int, (parent(xy) .- MD.o)./MD.d)

    # 边界检查保持不变
    if 0<i<MD.s[1] && 0<j<MD.s[2]
        # 直接使用 OffsetArray 的索引，不需要手动调整 offset
        n_vertices = @inbounds MD.T0[i, j, 0]::Int

        # 预获取底层数组以提高性能
        xy_data = parent(xy)::Vector{Float64}

        # 遍历顶点
        for vertex_idx in 1:n_vertices
            vertex = @inbounds MD.T0[i, j, vertex_idx]::Int
            n_tris = @inbounds MD.T1[6, vertex]::Int

            # 遍历三角形
            for tri_idx in 0:(n_tris-1)
                tri = @inbounds MD.T1[tri_idx, vertex]::Int

                # 使用视图避免分配
                tri_vertices = MD.T2[:, tri]
                tri_coords = MD.XY[:, tri_vertices]

                if is_inside(xy_data, parent(tri_coords))::Bool
                    return Vector{Int}(parent(tri_vertices))
                end
            end
        end
    end

    return Vector{Int}([-1,-1,-1])
end
function get_vertex_id_by_pos_new(xy::AbstractVector, MD::MetaData)::Vector{Int}
	i,j= floor.(Int, (parent(xy) .- MD.o)./MD.d)
	if 0<i<MD.s[1] && 0<j<MD.s[2]
		for vertex in MD.Tensor[1:MD.Tensor[0, j, i], j, i]
			for tri in MD.T1[0:MD.T1[6, vertex]-1, vertex]
				if is_inside(parent(xy), parent(MD.XY[:, MD.T2[:, tri]]))
					return parent(MD.T2[:, tri])
				end
			end
		end
	end
	error("出界!!!")
	return [-1,-1,-1]
end

function get_vertex_id_by_pos_type0(xy::AbstractVector, MD::MetaData)::Vector{Int}
    # 计算网格索引
    i, j = floor.(Int, (parent(xy) .- MD.o)./MD.d)

    # 边界检查
    if !(0 < i < MD.s[1] && 0 < j < MD.s[2])
        return Vector{Int}([-1, -1, -1])
    end

    # 明确类型标注
    n_vertices = @inbounds MD.T0[i, j, 0]::Int

    # 获取要用的底层数组以避免类型不稳定
    xy_data = parent(xy)::Vector{Float64}
    T2_data = parent(MD.T2)::Array{Int,2}
    XY_data = parent(MD.XY)::Array{Float64,2}

    # 遍历顶点
    for vertex_idx in 1:n_vertices
        vertex = @inbounds MD.T0[i, j, vertex_idx]::Int
        n_tris = @inbounds MD.T1[6, vertex]::Int

        # 遍历三角形
        for tri_idx in 0:(n_tris-1)
            tri = @inbounds MD.T1[tri_idx, vertex]::Int

            # 直接使用底层数组而不是视图
            # 获取三角形顶点编号
            tri_vertices = @view T2_data[:, tri]
            # 获取三角形顶点坐标
            tri_coords = @view XY_data[:, tri_vertices]

            if is_inside(xy_data, parent(tri_coords))
                # 直接使用类型转换确保返回类型正确
                return Vector{Int}(tri_vertices)
            end
        end
    end

    return Vector{Int}([-1, -1, -1])
end

"""
类型稳定
"""
function get_vertex_id_by_pos_type1(xy::AbstractVector, MD::MetaData)::Vector{Int}
    # 计算网格索引
    i, j = floor.(Int, (parent(xy) .- MD.o)./MD.d)

    # 边界检查
    if !(0 < i < MD.s[1] && 0 < j < MD.s[2])
        return Vector{Int}([-1, -1, -1])
    end

    # 预获取底层数组数据
    xy_data = parent(xy)::Vector{Float64}
    T0_data = parent(MD.T0)::Array{Int,3}  # 添加了 T0 的底层数组
    T1_data = parent(MD.T1)::Array{Int,2}  # 添加了 T1 的底层数组
    T2_data = parent(MD.T2)::Matrix{Int}
    XY_data = parent(MD.XY)::Matrix{Float64}

    # 获取 OffsetArray 的偏移
    T0_offsets = MD.T0.offsets
    T1_offsets = MD.T1.offsets

    # 计算实际索引
    i_idx = i - T0_offsets[1]
    j_idx = i - T0_offsets[2]

    # 获取顶点数量
    n_vertices = T0_data[i_idx, j_idx, 0 - T0_offsets[3]]::Int

    # 遍历顶点
    for vertex_idx in 1:n_vertices
        vertex = T0_data[i_idx, j_idx, vertex_idx - T0_offsets[3]]::Int
        n_tris = T1_data[6 - T1_offsets[1], vertex - T1_offsets[2]]::Int

        # 遍历三角形
        for tri_idx in 0:(n_tris-1)
            tri = MD.T1[tri_idx, vertex]::Int

            # 直接使用底层数组的视图
            tri_vertices = @view T2_data[:, tri]
            tri_coords = @view XY_data[:, tri_vertices]

            # 确保 is_inside 返回 Bool
            if is_inside(xy_data, parent(tri_coords))
                return Vector{Int}(tri_vertices)
            end
        end
    end

    return Vector{Int}([-1, -1, -1])
end
"""
类型稳定
"""
function get_vertex_id_by_pos_type(xy::AbstractVector, MD::MetaData)::Vector{Int}
    # 计算网格索引
    i, j = floor.(Int, (parent(xy) .- MD.o)./MD.d)

    # 边界检查
    if !(0 < i < MD.s[1] && 0 < j < MD.s[2])
        return Vector{Int}([-1, -1, -1])
    end

    # 预获取底层数组数据
    xy_data = parent(xy)::Vector{Float64}
    T0_data = parent(MD.T0)::Array{Int,3}  # 添加了 T0 的底层数组
    T1_data = parent(MD.T1)::Array{Int,2}  # 添加了 T1 的底层数组
    T2_data = parent(MD.T2)::Matrix{Int}
    XY_data = parent(MD.XY)::Matrix{Float64}

    # 获取 OffsetArray 的偏移
    T0_offsets = MD.T0.offsets
    T1_offsets = MD.T1.offsets

    # 计算实际索引
    i_idx = i - T0_offsets[1]
    j_idx = i - T0_offsets[2]

    # 获取顶点数量
    n_vertices = T0_data[i_idx, j_idx, 0 - T0_offsets[3]]::Int

    # 遍历顶点
    for vertex_idx in 1:n_vertices
        vertex = T0_data[i_idx, j_idx, vertex_idx - T0_offsets[3]]::Int
        n_tris = T1_data[6 - T1_offsets[1], vertex - T1_offsets[2]]::Int

        # 遍历三角形
        for tri_idx in 0:(n_tris-1)
            tri = T1_data[tri_idx - T1_offsets[1], vertex - T1_offsets[2]]::Int

            # 直接使用底层数组的视图
            tri_vertices = @view T2_data[:, tri]
            tri_coords = @view XY_data[:, tri_vertices]

            # 确保 is_inside 返回 Bool
            if is_inside(xy_data, parent(tri_coords))::Bool
                return Vector{Int}(tri_vertices)
            end
        end
    end

    return Vector{Int}([-1, -1, -1])
end
"""
通过点的坐标判断该点在磁面第几层，
层数从 1 开始数，即，最小层数为 1
"""
function get_psi_by_pos(xy::AbstractVector, MD::MetaData)::Int # 出 ITER 非截面
    tmp = get_vertex_id_by_pos(xy, MD)
    if tmp == [-1,-1,-1]
        return -1
    end
	return tmp .|> i2ψ |> maximum
end

function get_psi_by_pos(XY::Matrix, MD::MetaData) # 出 ITER 非截面
	ψs = zeros(size(XY, 2))
	for i in 1:size(XY, 2)
		ψs[i] = get_psi_by_pos(XY[:, i], MD::MetaData)
	end
	return ψs/100 # TODO 这里若 199 层，就 / 199
end
# %%

"""
Data.h5 大小: 多少 G
"""
data_size(p::Particles)=p.l*p.n*8*(11+3)/1024^3

"""
笛卡尔坐标转为极坐标
"""
function cart2polar(x::T,y::T) where T<:Real
	R = sqrt(x^2 + y^2);
	if R == 0
		ϕ = 0
		println("Singular point occors when CART to CYLD")
	end
	if y >= 0 
		ϕ = acos(x/R) 
	else
		ϕ = 2π-acos(x/R)
	end
	return [R, ϕ]
end

"""
笛卡尔坐标转为柱坐标
"""
function cart2cyld(x::T,y::T,z::T) where T<:Real
    R, ϕ=cart2polar(x,y)
	return [R, ϕ, z]
end
cart2cyld(v::Vector)=length(v) == 3 ? cart2cyld(v...) : error("dim != 3")

function cart2cyld(X::Matrix)
	rϕz = similar(X)
	if size(X,2) == 3
		for i = 1:size(X,1)
			rϕz[i, :] = cart2cyld(X[i, :])
		end
	else
		error("size 不为 3!!!")
	end
	return rϕz
end

function plot_tri(x::Vector,v::Matrix)
	plt.scatter(v[1, 1], v[2, 1], c="r")
	plt.scatter(v[1, 2], v[2, 2], c="g")
	plt.scatter(v[1, 3], v[2, 3], c="b")
	plt.plot(v[1, 1:2], v[2, 1:2])
	plt.plot(v[1, 2:3], v[2, 2:3])
	plt.plot(v[1, [3, 1]], v[2, [3,1]])
	plt.scatter(x[1], x[2], c="k")
	plt.show()
end


interpolateB(w::Vector{Float64}, subB::AbstractArray, id::Vector{Int})=sum(parent(subB[:, id]).*w', dims=2) |> vec

function cart2polar(p::Particle; R0=6.2)
    rϕ = zeros(size(p.X, 1), 2)
	for i in 1:size(p.X, 1)
        x,y,z=p.X[i, :]
        R = sqrt(x^2 + y^2) - R0
		rϕ[i, :] = cart2polar(R, z)
	end
	return rϕ
end
function cart2cyld(p::Particle)
	rϕz = similar(p.X)
	for i in 1:size(p.X, 1)
		rϕz[i, :] = cart2cyld(p.X[i, :])
	end
	return rϕz
end

"""
计算环向角
"""
function countϕ(ϕ::Vector{Float64})
	count = 0
	for i in 2:length(ϕ)
		if ϕ[i] < ϕ[i-1]
			count += 1
		end
	end
	return count
end
function countϕ(p::Particle)
	rϕz=cart2cyld(p)
	countϕ(rϕz[:, 2])
end

## %% gyro center {{{

#"""
#方均根
#"""
#rsq(B::AbstractMatrix) = sum(B.^2, dims=2) .|> sqrt |> vec

#"""
#用回旋半径
#计算导心
#``R = mv/qB``
#``R = mv \\times B / (qB^2)``
#	R = p × B / (u.q * sum(B.^2))
#"""
#function gyro_R(p::Particle)::Matrix
#	R = similar(p.P)
#	for i = 1:size(R, 1)
#		P = p.P[i, :]
#		B = p.B[i, :]
#		R[i, :] = P × B / (u.q * sum(B.^2))
#	end
#	return R
#end

#function center(p::Particle)::Matrix
#	p.X + gyro_R(p)
#end


"""
是否是通行粒子
"""
is_passage(p::Particle)::Bool = p|> p->eachrow(p.B) .⋅ eachrow(p.P) .|> sign |> unique |> length == 1

##}}}


end
