#! /usr/bin/env -S julia --color=yes --startup-file=no
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#=
    test
    Copyright © 2024 <junyixu0@gmail.com>

    Distributed under terms of the MIT license.
=#
# %%
const Δt = 1.0e-3
const N = 2 # number of particles
const TotalSteps = 1000 # total steps
# %%

include("PtcStruct.jl")
using .PtcStruct
include("Pushers.jl")
using .Pushers:pusher_gravity as pusher
using PyCall
# %%
@pyimport matplotlib.pyplot as plt
@pyimport celluloid
# %%
function init!(ps::Vector{Particle})
    R = 1
    x = [R 0]
    v = [0 0.5/sqrt(R)]
    m = 1.0
    ps[1] = Particle([x ;zeros(TotalSteps, 2)], [v ;zeros(TotalSteps, 2)], m)
    ps[2] = Particle([-x ;zeros(TotalSteps, 2)], [-v ;zeros(TotalSteps, 2)], m)
end
# %%

function push_ptcs!(ps::Vector{Particle})
    for s = 1:TotalSteps 
        for i = 1:N
            ps[i].X[s+1, :], ps[i].V[s+1, :] = pusher(ps, i, s)
        end
    end
end

function anim(ps::Vector{Particle})
    fig = plt.figure(dpi = 300)
    camera = celluloid.Camera(fig)

    for s in 1:10:TotalSteps
        # for i = 1:N
        #     plt.scatter(ps[i].X[s, 1], ps[i].X[s, 2])
        # end
        plt.scatter(ps[1].X[s, 1], ps[1].X[s, 2], c="r")
        plt.plot(ps[1].X[1:s, 1], ps[1].X[1:s, 2], "r--")
        plt.scatter(ps[2].X[s, 1], ps[2].X[s, 2], c="b")
        plt.plot(ps[2].X[1:s, 1], ps[2].X[1:s, 2], "b--")
        plt.axis("equal")
        camera.snap()
    end

    animation = camera.animate(blit=false, interval=100, repeat=false)
    animation.save("nbodies.mp4",
                   dpi=300,
                   savefig_kwargs=Dict("pad_inches"=>"tight"))
    # plt.show()

end
# %%

function main()
    ps=Vector{Particle}(undef, N) # Particles
    init!(ps)
    @time push_ptcs!(ps)
    @time anim(ps)
end
# %%
main()
