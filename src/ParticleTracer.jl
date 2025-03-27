"""
    ParticleTracer

高能粒子追踪器
"""
module ParticleTracer
using Distributed


include("Constants.jl")
include("UserInputs.jl")
include("PtcStruct.jl")
include("Fields.jl")
include("Pushers.jl")
include("DataIO.jl")


using .PtcStruct
using .PtcStruct:MagneticParticleData, EMParticleData
pusher = @eval Pushers.$(UserInputs.pusher)
get_fields = @eval Fields.$(UserInputs.field)
using .UserInputs: TotalSteps, SavePerNSteps, is_data_saving_on, is_merge_process_files_on
using .Constants
using HDF5
using LinearAlgebra: ⋅, norm
include("initialization.jl")
using .DataIO

# %%
function init_ptc_data(x0::AbstractVector, p0::AbstractVector, N::Int)
    x0 = reshape(x0, 3, 1)
    p0 = reshape(p0, 3, 1)
    X = [x0 zeros(3, N-1)]
    P = [p0 zeros(3, N-1)]
    B = zeros(3, N)
    
    if UserInputs.use_electric_field
        E = zeros(3, N)
        return ParticleData(X, P, B, E)
    else
        return ParticleData(X, P, B)
    end
end
function save(ptc_data::MagneticParticleData, ptc::MagneticParticle, i::Int)
    ptc_data.X[:, i] .= ptc.X
    ptc_data.P[:, i] .= ptc.P
    ptc_data.B[:, i] .= ptc.B
end
function save(ptc_data::EMParticleData, ptc::EMParticle, i::Int)
    ptc_data.X[:, i] .= ptc.X
    ptc_data.P[:, i] .= ptc.P
    ptc_data.B[:, i] .= ptc.B
    ptc_data.E[:, i] .= ptc.E
end
function push_ptc!(ptc::MagneticParticle)
    B, E = get_fields(ptc.X...)
    
    # 检查是否出界
    if all(iszero, B) && all(iszero, E)
        return false
    end
    
    xx, pp = pusher(ptc.X, ptc.P, E, B)
    
    for i in 1:3
        ptc.X[i] = xx[i]
        ptc.P[i] = pp[i]
        ptc.B[i] = B[i]
    end
    return true
end
function push_ptc!(ptc::EMParticle)
    B, E = get_fields(ptc.X...)
    # 检查是否出界
    if all(iszero, B) && all(iszero, E)
        return false
    end
    xx, pp = pusher(ptc.X, ptc.P, E, B)
    
    # 逐元素更新，而不是整体赋值
    for i in 1:3
        ptc.X[i] = xx[i]
        ptc.P[i] = pp[i]
        ptc.B[i] = B[i]
        ptc.E[i] = E[i]
    end
    return true
end

# %%
# function anim(ps::Vector{Particle})
#     fig = plt.figure(dpi = 300)
#     camera = celluloid.Camera(fig)

#     for s in 1:10:TotalSteps
#         # for i = 1:N
#         #     plt.scatter(ps[i].X[s, 1], ps[i].X[s, 2])
#         # end
#         plt.scatter(ps[1].X[s, 1], ps[1].X[s, 2], c="r")
#         plt.plot(ps[1].X[1:s, 1], ps[1].X[1:s, 2], "r--")
#         plt.scatter(ps[2].X[s, 1], ps[2].X[s, 2], c="b")
#         plt.plot(ps[2].X[1:s, 1], ps[2].X[1:s, 2], "b--")
#         plt.axis("equal")
#         camera.snap()
#     end

#     animation = camera.animate(blit=false, interval=100, repeat=false)
#     animation.save("nbodies.mp4",
#                    dpi=300,
#                    savefig_kwargs=Dict("pad_inches"=>"tight"))
#     # plt.show()

# end

# %%
"""
根据配置生成初始位置和动量
"""
function generate_initial_conditions(init_type::Symbol)
    if init_type == :single
        # 单粒子，使用确定的初始值
        return [UserInputs.x0], [UserInputs.p0]
    elseif init_type == :parabolic_torus
        # 从环面抛物线分布采样多粒子
        x0_list = [SetParticlePosition_ParabolicTorus(UserInputs.r_max) for _ in 1:UserInputs.N]
        
        # 对每个位置计算对应的动量
        p0_list = [SetParticleMomentum_Gyrocenter(x..., get_fields(x...)[1]) 
                   for x in x0_list] # p = f(x,y,z,Bx,By,Bz), 其中 Bx,By,Bz 是磁场, 无电场
        
        return x0_list, p0_list
    else
        error("未知的初始化类型: $init_type")
    end
end

function simulate_particle!(ptc_data, ptc, n, save_per_n_steps=SavePerNSteps)
    num_outputs = 10  # 期望的输出次数
    output_interval = (TotalSteps - 1) ÷ num_outputs
    # 预先计算输出步数
    output_steps = Set(output_interval:output_interval:TotalSteps-1)
    
    for i in 1:TotalSteps-1
        if !(push_ptc!(ptc))
            return i  # 返回出界时的步数
        end

        if myid() == 1 && n==1 && i in output_steps
            progress = round(i / (TotalSteps-1) * 100, digits=1)
            println("模拟进度: $(progress)% (步骤 $i / $(TotalSteps))")
        end
        
        if iszero(i % save_per_n_steps) && i != TotalSteps
            save(ptc_data, ptc, i ÷ save_per_n_steps + 1)
        end
    end
    return nothing  # 粒子未出界
end

function initialize_particle(x0::Vector{T}, p0::Vector{T}) where T<:AbstractFloat
    B0, E0 = get_fields(x0...)
    UserInputs.use_electric_field ? Particle(x0, p0, B0, E0) : Particle(x0, p0, B0)
end

function main()
    t_start = time()
    t_io = 0.0  # IO操作累计时间
    
    # 根据进程数调整批处理大小
    num_workers = nworkers()
    optimal_batch_size = max(1000 ÷ num_workers, 100)  # 根据进程数调整批大小
    save_config = SaveConfig(batch_size=optimal_batch_size)
    data_length = TotalSteps ÷ SavePerNSteps
    
    # 生成并处理所有粒子
    # x0_list, p0_list = generate_initial_conditions(UserInputs.init_type)
    
    # 创建进度跟踪变量
    last_sync_step = 0
    # 使用全局共享的RemoteChannel
    
    @sync @distributed for n in 1:UserInputs.N
        x0 = SetParticlePosition_ParabolicTorus(UserInputs.r_max)
        p0 = SetParticleMomentum_Gyrocenter(x0..., get_fields(x0...)[1])
        
        ptc = initialize_particle(x0, p0)
        ptc_data = init_ptc_data(x0, p0, data_length)
        
        if (escape_step = simulate_particle!(ptc_data, ptc, n)) !== nothing
            @info "粒子 $n 在第 $escape_step 步出界"
        end
        # t_sim = time() - t_start
        # println("进程 myid():\t模拟时间: $(round(t_sim, digits=4)) 秒")
        save_particle_data(ptc_data, n, x0, p0, get_fields(x0...)[1], save_config)
        # println("进程 myid():\t写入时间: $(round((time() - t_start), digits=4)) 秒")
    end
    t_sim = time() - t_start
    # 收集所有出界粒子信息
    n_total = UserInputs.N

    # 记录最终IO操作时间
    t_io_start = time()
    if is_data_saving_on && is_merge_process_files_on
        merge_process_files(save_config)
        create_index_file(save_config)
    end
    t_io += time() - t_io_start

    t_total = time() - t_start
    println("进程 $(myid()):")
    println("\t总计算时间: $(round(t_total, digits=4)) 秒")
    println("\t模拟时间: $(round(t_sim, digits=4)) 秒 ($(round(t_sim/t_total*100, digits=4))%)")
    println("\tIO时间: $(round(t_io, digits=4)) 秒 ($(round(t_io/t_total*100, digits=4))%)")
end
# %%

#=
f = h5open("ptc_data4.h5", "r")
X = read(f, "X")

    x = X[1, :]
    y = X[2, :]
    z = X[3, :]
    x = x * u.x
    y = y * u.x
    z = z * u.x
    R = sqrt.(x.^2 + y.^2)


    x = ptc_data.X[1, :]
    y = ptc_data.X[2, :]
    z = ptc_data.X[3, :]
    x = x * u.x
    y = y * u.x
    z = z * u.x
    R = sqrt.(x.^2 + y.^2)

    plt.plot(R, z)
    # plt.axis("equal")
    plt.show()
=#

    # ax=plt.axes(projection="3d")
    # ax.plot3D(x, y, z)
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # ax.set_zlabel("Z")
    # plt.show()


end # module ParticleTracer
