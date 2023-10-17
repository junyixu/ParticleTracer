module Plots
using HDF5
using PyCall
using LaTeXStrings
export plot_phase, plot_E

@pyimport matplotlib.pyplot as plt


function plot_phase(filename::String; color::String="b")
	h5open(filename, "r") do f
		x = f["x"] |> read
		p = f["p"] |> read
		plt.scatter(x,p, c=color, s=1)
	end
end

function plot_E(filename::String; color::String="b")
	h5open(filename, "r") do ff
		x = ff["x"] |> read
		p = ff["p"] |> read
		plt.plot(0.5 .* p.^2 + 0.5 .* x.^2, color=color, linewidth=1)
	end
end

end
