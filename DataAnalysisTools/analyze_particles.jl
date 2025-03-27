using HDF5
using Dates
using Plots
using LinearAlgebra
using Statistics
include("../src/PtcStruct.jl")
using .PtcStruct
include("../src/Constants.jl")
using .Constants


# %%
"""
分析结果的数据结构
"""
struct TrajectoryAnalysis
    R::Vector{Float64}        # 径向位置
    Z::Vector{Float64}        # 纵向位置
    energy::Vector{Float64}   # 能量
    raw_data::AbstractParticleData
end

"""
粒子集合容器
"""
struct ParticleEnsemble
    timestamp::String
    particles::Vector{AbstractParticleData}
    metadata::Dict{String, Any}  # 存储批次级别的元数据
end

# 添加索引操作符支持
Base.getindex(pe::ParticleEnsemble, i::Int) = pe.particles[i]
Base.length(pe::ParticleEnsemble) = length(pe.particles)
Base.iterate(pe::ParticleEnsemble) = iterate(pe.particles)
Base.iterate(pe::ParticleEnsemble, state) = iterate(pe.particles, state)
Base.length(ptcl::AbstractParticleData) = length(ptcl.X[1, :])

"""
查找最新的索引文件并返回其时间戳
"""
function get_latest_timestamp()
    # 筛选出索引文件
    index_files = filter(x -> startswith(x, "index_"), readdir())
    if isempty(index_files)
        error("未找到索引文件")
    end
    
    # 从文件名提取时间戳
    timestamps = map(f -> match(r"index_(\d+_\d+).h5", f).captures[1], index_files)
    
    # 返回最新的时间戳
    return sort(timestamps)[end]
end

"""
从指定时间戳的批次文件中读取粒子数据
"""
function read_particle_data(timestamp::AbstractString, particle_id::Int)::AbstractParticleData
    # 确定批次号
    batch_number = ceil(Int, particle_id/1000)
    filename = "particles_$(timestamp)_batch$(batch_number).h5"
    
    h5open(filename, "r") do f
        particle = f["particle_$(mod1(particle_id, 1000))"]
        traj = particle["trajectory"]
        
        pos = read(traj["position"])
        mom = read(traj["momentum"])
        B = read(traj["magnetic_field"])
        
        if haskey(traj, "electric_field")
            E = read(traj["electric_field"])
            return EMParticleData(pos, mom, B, E)
        else
            return MagneticParticleData(pos, mom, B)
        end
    end
end

# %%
"""
A wrapper struct that combines particle data with its unit system.
This allows for:
1. Encapsulation of both data and its physical interpretation
2. Cleaner function interfaces
3. Consistent unit handling
while maintaining the flexibility to use different unit systems when needed.
"""
struct UnitParticleData <: AbstractParticleData
    data::AbstractParticleData    # Raw particle data (positions, momenta, fields)
    units::Unit                   # Physical unit system for interpretation
end

# Delegate common fields and methods to the underlying particle data
Base.length(upd::UnitParticleData) = length(upd.data)
Base.getproperty(upd::UnitParticleData, sym::Symbol) = 
    sym ∈ (:data, :units) ? getfield(upd, sym) : getproperty(upd.data, sym)
# %%

# 定义简单的Trait
abstract type HasUnits end
struct WithUnits <: HasUnits end
struct NoUnits <: HasUnits end

# Trait接口
has_units(::AbstractParticleData) = NoUnits()
has_units(::UnitParticleData) = WithUnits()
get_units(ptc::UnitParticleData) = ptc.units

# 复用原有的calculate_energy，但增加Trait分派
function calculate_energy(ptcl::AbstractParticleData; i = length(ptcl))
    calculate_energy(ptcl, has_units(ptcl), i)
end

function calculate_energy(ptcl::AbstractParticleData, ::NoUnits, i)
    error("需要带单位的粒子数据")
end

function calculate_energy(ptcl::AbstractParticleData, ::WithUnits, i)
    calculate_energy(ptcl.data, get_units(ptcl), i)
end

# %%

# Example usage:
# Single particle processing
particle_data=ps[1]
particle_with_units = UnitParticleData(particle_data, u)
energy = calculate_energy(particle_with_units)

particle_data=ps[1]
particle_with_units = UnitParticleData(particle_data, u)
particle_with_units.t
energy = calculate_energy(particle_with_units)

# Batch processing with shared units
ensemble_with_units = map(p -> UnitParticleData(p, u), particles)
energies = calculate_energy.(ensemble_with_units)

"""
将笛卡尔坐标转换为柱坐标(R,φ,Z)
"""
function cart2cyl(x::Vector{Float64})
    R = sqrt(x[1]^2 + x[2]^2)
    φ = atan(x[2], x[1])
    return [R, φ, x[3]]
end
# %%

function cart2cyl(ptcl::AbstractParticleData)
    cart2cyl(ptcl, has_units(ptcl))
end

function cart2cyl(pos::Matrix)
    R = sqrt.(pos[1,:].^2 + pos[2,:].^2)
    φ = atan.(pos[2,:], pos[1,:])
    return (R, φ, pos[3,:])
end

function cart2cyl(ptcl::AbstractParticleData, ::NoUnits)
    R, φ, z = cart2cyl(ptcl.X)
    return R, φ, z
end

function cart2cyl(ptcl::AbstractParticleData, ::WithUnits)
    R, φ, z = cart2cyl(ptcl.X)
    R *= get_units(ptcl).x
    z *= get_units(ptcl).x
    return R, φ, z
end

# %%

p2 = read_particle_data(timestamp, 2)

R, φ, z=cart2cyl(ps[1])
plot(R*u.x,z*u.x)

# %%

pos = ptc_data.X * u.x
mom = ptc_data.P * u.p
cyl_coords = [cart2cyl(pos[:,i]) for i in 1:size(pos,2)]
R = getindex.(cyl_coords, 1)
Z = getindex.(cyl_coords, 3)

# %%
"""
分析单个粒子的轨迹
"""
function analyze_trajectory(ptc_data::AbstractParticleData, u)
    # 转换为物理单位
    pos = ptc_data.X * u.x
    mom = ptc_data.P * u.p
    
    # 计算能量
    energy = [calculate_energy(mom[:,i], u) for i in 1:size(mom,2)]
    
    # 转换为柱坐标
    cyl_coords = [cart2cyl(pos[:,i]) for i in 1:size(pos,2)]
    R = getindex.(cyl_coords, 1)
    Z = getindex.(cyl_coords, 3)
    
    return TrajectoryAnalysis(R, Z, energy, ptc_data)
end

# %%
"""
绘制粒子轨迹
"""
function plot_trajectory(analysis::TrajectoryAnalysis)
    # R-Z平面轨迹图
    p1 = plot(analysis.R, analysis.Z,
        xlabel="R (m)",
        ylabel="Z (m)",
        title="Particle Trajectory",
        aspect_ratio=:equal,
        label="trajectory"
    )
    
    # 能量随时间变化
    p2 = plot(analysis.energy,
        xlabel="Time Step",
        ylabel="Energy (keV)",
        title="Particle Energy",
        label="energy"
    )
    
    # 组合图
    plot(p1, p2, layout=(2,1), size=(800,1000))
end

"""
从指定时间戳读取所有粒子数据
"""
function load_particle_ensemble(timestamp::AbstractString; max_particles::Int=0)
    # 读取索引文件获取元数据
    index_filename = "index_$timestamp.h5"
    metadata = h5open(index_filename, "r") do f
        Dict(String(name) => read_attribute(f, name) 
            for name in keys(attrs(f)))
    end
    
    total_particles = metadata["total_particles"]
    batch_size = metadata["batch_size"]
    
    # 如果指定了最大粒子数，则限制读取数量
    n_particles = max_particles > 0 ? min(max_particles, total_particles) : total_particles

    # 预分配粒子数组
    particles = Vector{AbstractParticleData}(undef, n_particles)
    
    # 按批次读取粒子数据
    for particle_id in 1:n_particles
        particles[particle_id] = read_particle_data(timestamp, particle_id)
        
        # 显示进度
        if particle_id % 100 == 0
            println("已读取 $particle_id / $n_particles 个粒子")
        end
    end
    
    return ParticleEnsemble(timestamp, particles, metadata)
end

"""
分析粒子集合
"""
struct EnsembleAnalysis
    trajectories::Vector{TrajectoryAnalysis}
    statistics::Dict{String, Any}  # 存储统计信息
end

function analyze_ensemble(ensemble::ParticleEnsemble, u)
    # 分析每个粒子的轨迹
    trajectories = [analyze_trajectory(ptc, u) for ptc in ensemble.particles]
    
    # 计算统计量
    statistics = Dict{String, Any}()
    
    # 计算平均能量随时间的变化
    mean_energy = mean([t.energy for t in trajectories])
    std_energy = std([t.energy for t in trajectories])
    
    # 计算径向分布
    final_R = [t.R[end] for t in trajectories]
    R_histogram = fit(Histogram, final_R)
    
    statistics["mean_energy"] = mean_energy
    statistics["std_energy"] = std_energy
    statistics["R_histogram"] = R_histogram
    
    return EnsembleAnalysis(trajectories, statistics)
end

"""
绘制集合分析结果
"""
function plot_ensemble_analysis(analysis::EnsembleAnalysis)
    # 轨迹叠加图
    p1 = plot(title="Particle Trajectories (R-Z)")
    for traj in analysis.trajectories
        plot!(p1, traj.R, traj.Z, 
            alpha=0.1,  # 透明度
            label="",   # 不显示图例
            xlabel="R (m)",
            ylabel="Z (m)",
            aspect_ratio=:equal
        )
    end
    
    # 平均能量演化
    p2 = plot(
        analysis.statistics["mean_energy"],
        ribbon=analysis.statistics["std_energy"],
        title="Mean Particle Energy",
        xlabel="Time Step",
        ylabel="Energy (keV)",
        label="Mean ± Std"
    )
    
    # 最终径向分布
    hist = analysis.statistics["R_histogram"]
    p3 = plot(hist, 
        title="Final Radial Distribution",
        xlabel="R (m)",
        ylabel="Count",
        label="Particles"
    )
    
    # 组合图
    plot(p1, p2, p3, layout=(2,2), size=(1000,1000))
end

# 主程序示例
function main()
    timestamp = get_latest_timestamp()
    println("分析时间戳: $timestamp 的数据")
    
    # 读取粒子集合（这里限制为前1000个粒子作为示例）
    ensemble = load_particle_ensemble(timestamp, max_particles=1000)
    println("已加载 $(length(ensemble.particles)) 个粒子")
    
    # 设置物理单位
    u = Unit(5.18, "alpha")
    
    # 分析粒子集合
    analysis = analyze_ensemble(ensemble, u)
    
    # 绘制并保存分析结果
    p = plot_ensemble_analysis(analysis)
    savefig(p, "ensemble_analysis_$(timestamp).pdf")
    
    println("分析完成，图像已保存")
end

# 运行主程序
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end 
