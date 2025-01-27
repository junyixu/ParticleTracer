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


@everywhere begin
include("Fields.jl")
include("PtcStruct.jl")
include("UserInputs.jl")
include("Pushers.jl")
include("Constants.jl")
include("DataIO.jl")
end


@everywhere begin
using .PtcStruct
pusher = @eval Pushers.$(UserInputs.pusher)
using .UserInputs: TotalSteps, SavePerNSteps
using .Constants
using HDF5
using LinearAlgebra: ⋅, norm
include("initialization.jl")
end

@everywhere using .DataIO

# %%
@everywhere begin
function init_ptc_data(x0::AbstractVector, p0::AbstractVector, N::Int)
    x0 = reshape(x0, 3, 1)
    p0 = reshape(p0, 3, 1)
    X = [x0 zeros(3, N-1)]
    P = [p0 zeros(3, N-1)]
    B = zeros(3, N)
    ptc_data = ParticleData(X, P, B)
    return ptc_data
end
function save(ptc_data::ParticleData, ptc::Particle, i::Int)
    ptc_data.X[:, i] .= ptc.X
    ptc_data.P[:, i] .= ptc.P
    ptc_data.B[:, i] .= ptc.B
end
function push_ptc!(ptc)
    x = ptc.X # vector x
    p = ptc.P # vector p
    B = Fields.tokamak(x..., 1.0) # q = 1.0
    xx, pp = pusher(x, p, [0.0, 0, 0], B)
    ptc.X .= xx
    ptc.P .= pp
    ptc.B .= B
end
end # @everywhere

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
    # Create save configuration
    save_config = SaveConfig(batch_size=1000)
    
    @sync @distributed for n = 1:UserInputs.N
        # Initialize particle parameters
        x0 = SetParticlePosition_ParabolicTorus(UserInputs.a)
        B0 = Fields.tokamak(x0..., 1.0)
        p0 = SetParticleMomentum_Gyrocenter(x0..., B0) 
        ptc = Particle(x0, p0, B0)
        data_length = Int(TotalSteps/SavePerNSteps)
        ptc_data = init_ptc_data(x0, p0, data_length)

        # Main computation loop
        for i in 2:TotalSteps-1
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
