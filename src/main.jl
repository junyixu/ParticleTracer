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
include("Pushers.jl")
include("UserInputs.jl")
using .UserInputs: TotalSteps
include("Constants.jl")
using .Constants

# %%
function init_ptc(x0::AbstractVector, p0::AbstractVector)
    x0 = reshape(x0, 1, 3)
    p0 = reshape(p0, 1, 3)
    X = [x0; zeros(TotalSteps, 3)]
    P = [p0; zeros(TotalSteps, 3)]
    B = zeros(TotalSteps+1, 3)
    ptc = Particle(X, P, B)
    return ptc
end

function push_ptc!(ptc)
    for s in 1:TotalSteps
        x = ptc.X[s, :] # vector x
        p = ptc.P[s, :] # vector p
        B = Fields.tokamak(x..., 1.0) # q = 1.0
        ptc.B[s+1, :] = B
        ptc.X[s+1, :], ptc.P[s+1, :] = Pushers.boris(x, p, [0, 0, 0], B)
    end
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
    ptc = init_ptc(Constants.x0,Constants.p0)
    push_ptc!(ptc)

    x = ptc.X[:, 1]
    y = ptc.X[:, 2]
    z = ptc.X[:, 3]
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
