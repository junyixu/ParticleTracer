# %%
include("Pushers.jl")
using .Pushers
include("Param.jl")
using .Param
include("Plots.jl")
using .Plots
using HDF5
using PyCall
using LaTeXStrings
# %%
const output_dir = "output/"

# %%
@pyimport matplotlib.pyplot as plt

# %%
function main(P::Parameters)
	x=zeros(P.nsteps)
	p=zeros(P.nsteps)

	x[1] = P.x0
	p[1] = P.p0

	for j = 1:P.nfiles
		for i = 2:P.nsteps
			x[i],p[i]=P.pusher(x[i-1],p[i-1])
		end
		x[1],p[1] = P.pusher(x[end], p[end])
		h5open("$(output_dir*P.filename)$j.h5", "w") do f
			f["x"]=x
			f["p"]=p
		end
	end
end
# %%
p=Param.p
main(p)

# %%

plot_phase("$(output_dir*p.filename)1.h5", color="b")
plot_phase("$(output_dir*p.filename)$(p.nfiles).h5", color="r")
plt.axis("equal")
plt.savefig("./figures/$(p.filename).pdf", bbox_inches="tight")
plt.show()
plot_E("$(output_dir*p.filename)1.h5", color="b")
plot_E("$(output_dir*p.filename)$(p.nfiles).h5", color="r")
plt.savefig("./figures/$(p.filename)E.pdf", bbox_inches="tight")
plt.show()
