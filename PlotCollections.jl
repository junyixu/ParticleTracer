module PlotCollections
using PyPlot

function p2scatter(p)
    x,y=p
    scatter(x,y, marker="p")
end

end # module
