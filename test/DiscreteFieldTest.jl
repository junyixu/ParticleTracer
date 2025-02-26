module DiscreteFieldTest
using Test

include("../src/DiscreteFields.jl")
include("../src/Fields.jl")
import .DiscreteFields:discrete

include("../PlotCollections.jl")
import .PlotCollections: p2scatter, savefig, scatter
include("/home/junyi/.julia/dev/Ptcs/src/Ptcs.jl")
import .Ptcs: get_vertex_id_by_pos!, MetaData, cart2cyld

# %%

const RESOURCE="/home/junyi/WorkSpace/apt_mhd/four_plus_100_ITER_100/"
MD = MetaData(RESOURCE*"tearing_mode_3tables_compress.h5")

"""
用三角形顶点编号 (vertex) 画一个平面上的点
"""
v2scatter(v::Int) = scatter(MD.XY[:, v]...)

x,y,z=(1.57941 , 5.24294 , -1.29892)

@time discrete(x,y,z)
@time discrete(x,y,z)

MD.T0 |> typeof # The type has been confirmed to be complete
# TODO 这样不方便，我还是想用昨天在 b 站上看到的用 AI 补全
three=[19964, 20063,20163]
three=zeros(Int, 3)
𝐱=[1.0355045208512683, -1.7406639441061689]
𝐱=[0.08915424, 0.7061096]
if -1 == @time get_vertex_id_by_pos!(three,𝐱, MD)
    !isdir("./figures") && mkpath("./figures")
    three .|> v2scatter
    𝐱 |> p2scatter
    savefig("./figures/pentagon_in_triangle.png", bbox_inches="tight")
end
end # module
