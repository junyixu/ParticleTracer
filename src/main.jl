#! /usr/bin/env -S julia --color=yes --startup-file=no
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#=
    test particles
    Copyright © 2024 Junyi Xu <junyixu0@gmail.com>

    Distributed under terms of the MIT license.
=#


using MyPlots
using Debugger

include("Fields.jl")
include("PtcStruct.jl")
using .PtcStruct
include("UserInputs.jl")
include("Pushers.jl")
pusher = @eval Pushers.$(UserInputs.pusher)
using .UserInputs: TotalSteps, SavePerNSteps
include("Constants.jl")
using .Constants
using HDF5

# %%
function init_ptc_data(x0::AbstractVector, p0::AbstractVector, N::Int)
    x0 = reshape(x0, 3, 1)
    p0 = reshape(p0, 3, 1)
    X = [x0 zeros(3, N-1)]
    P = [p0 zeros(3, N-1)]
    B = zeros(3, N)
    ptc_data = ParticleData(X, P, B)
    return ptc_data
end
# %%
function save(ptc_data::ParticleData, ptc::Particle, i::Int)
    ptc_data.X[:, i] .= ptc.X
    ptc_data.P[:, i] .= ptc.P
    ptc_data.B[:, i] .= ptc.B
end
# %%
function push_ptc!(ptc)
    x = ptc.X # vector x
    p = ptc.P # vector p
    B = Fields.tokamak(x..., 1.0) # q = 1.0
    xx, pp = pusher(x, p, [0.0, 0, 0], B)
    ptc.X .= xx
    ptc.P .= pp
    ptc.B .= B
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

# %%
function main()

    ptc = Particle(Constants.x0, Constants.p0, [0.0, 0, 0])
    data_length = Int(TotalSteps/SavePerNSteps)
    ptc_data = init_ptc_data(Constants.x0,Constants.p0, data_length)

    for i in 2:TotalSteps
        push_ptc!(ptc)
        # i % SavePerNSteps == 0 && @bp
        i % SavePerNSteps == 0 && i != TotalSteps && save(ptc_data, ptc, Int(i/SavePerNSteps)+1)
    end

    h5open("ptc_data.h5", "w") do file
        write(file, "X", ptc_data.X)
        write(file, "P", ptc_data.P)
        write(file, "B", ptc_data.B)
    end

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

    # ax=plt.axes(projection="3d")
    # ax.plot3D(x, y, z)
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # ax.set_zlabel("Z")
    # plt.show()

end
# %%
main()
