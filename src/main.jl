#! /usr/bin/env -S julia --color=yes --startup-file=no
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#=
    test particles
    Copyright © 2024 Junyi Xu <junyixu0@gmail.com>

    Distributed under terms of the MIT license.
=#

# using MyPlots
# using Debugger
using Dates

using Distributed


# @everywhere begin
include("Constants.jl")
include("UserInputs.jl")
include("PtcStruct.jl")
include("Fields.jl")
include("Pushers.jl")
include("DataIO.jl")
# end


# @everywhere begin
using .PtcStruct
using .PtcStruct:MagneticParticleData, EMParticleData
pusher = @eval Pushers.$(UserInputs.pusher)
get_fields = @eval Fields.$(UserInputs.field)
using .UserInputs: TotalSteps, SavePerNSteps
using .Constants
using HDF5
using LinearAlgebra: ⋅, norm
include("initialization.jl")
# end

using .DataIO

# %%
# @everywhere begin
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
    xx, pp = pusher(ptc.X, ptc.P, E, B)
    
    # 逐元素更新，而不是整体赋值
    for i in 1:3
        ptc.X[i] = xx[i]
        ptc.P[i] = pp[i]
        ptc.B[i] = B[i]
    end
end
function push_ptc!(ptc::EMParticle)
    B, E = get_fields(ptc.X...)
    xx, pp = pusher(ptc.X, ptc.P, E, B)
    
    # 逐元素更新，而不是整体赋值
    for i in 1:3
        ptc.X[i] = xx[i]
        ptc.P[i] = pp[i]
        ptc.B[i] = B[i]
        ptc.E[i] = E[i]
    end
end
# end # @everywhere

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
Main simulation program for particle trajectory calculation.
Handles parallel computation and data saving.
"""
function main()
    start_time = time()  # 记录开始时间
    
    # Create save configuration
    save_config = SaveConfig(batch_size=1000)
    
    for n = 1:UserInputs.N
        # Initialize particle parameters
        x0 = UserInputs.x0
        B0,E0 = get_fields(x0...)
        p0 = UserInputs.p0
        if UserInputs.use_electric_field
            ptc = Particle(x0, p0, B0, E0)
        else
            ptc = Particle(x0, p0, B0)
        end
        data_length = Int(TotalSteps/SavePerNSteps)
        ptc_data = init_ptc_data(x0, p0, data_length)

        # Main computation loop
        for i in 1:TotalSteps-1
            if myid() == 1 && (i % round(Int, UserInputs.TotalSteps/10) == 0 || i == UserInputs.TotalSteps-1)
                println("##\tstep = $i")
            end
            push_ptc!(ptc)
            # Save intermediate results
            i % SavePerNSteps == 0 && i != TotalSteps && 
                save(ptc_data, ptc, Int(i/SavePerNSteps)+1)
        end

        # Each process saves to its temporary file
        save_particle_data(ptc_data, n, x0, p0, B0, save_config)
    end
    
    # After all computations, merge process files
    merge_process_files(save_config)
    # Create index file for the dataset
    create_index_file(save_config)

    end_time = time()  # 记录结束时间
    if myid() == 1  # 只在主进程上打印
        println("进程 $(myid()), 总计算时间: $(end_time - start_time) 秒")
    end
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


# %%
main()
