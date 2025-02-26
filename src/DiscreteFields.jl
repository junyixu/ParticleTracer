module DiscreteFields
export discrete

include("UserInputs.jl")
R0 = UserInputs.R0
function cycle(i::Integer, L::Integer)
    if i < 0
        return i + L
    elseif i >= L
        return i - L
    else
        return i
    end
end

using HDF5
using OffsetArrays

include("/home/junyi/.julia/dev/Ptcs/src/Ptcs.jl")
import .Ptcs: get_vertex_id_by_pos!, MetaData, cart2cyld
using ..UserInputs: use_electric_field

const Δϕ = 2π/18
const RESOURCE="/home/junyi/WorkSpace/apt_mhd/four_plus_100_ITER_100/"
MD = MetaData(RESOURCE*"tearing_mode_3tables_compress.h5")

function h5load(filename::String, obj::String)
    fid=h5open(filename, "r")
    data=read(fid[obj])
    close(fid)
    return data
end

B=h5load(RESOURCE*"BX_BY_BZ_0.h5", "B") # TODO `reshape(:, :, 3)` in julia should be more intuitive than offset in c language, don't you think?
B = OffsetArray(reshape(B, 3, :, 18), 0, -1, -1)

E=h5load(RESOURCE*"BX_BY_BZ_0.h5", "E")
E = OffsetArray(reshape(E, 3, :, 18), 0, -1, -1)

function interpolate𝐁(three::Vector{Int}, 𝐱, 𝐁_data)
    𝐯 = [parent(MD.XY[:, v]) for v in three]
    𝐁s = [𝐁_data[:, v] for v in three]

    # begin 计算重心坐标权重 𝐰
    den = ((𝐯[2][2] - 𝐯[3][2]) * (𝐯[1][1] - 𝐯[3][1]) +
           (𝐯[3][1] - 𝐯[2][1]) * (𝐯[1][2] - 𝐯[3][2]))

    𝐰 = zeros(3)
    𝐰[1] = ((𝐯[2][2] - 𝐯[3][2]) * (𝐱[1] - 𝐯[3][1]) +
             (𝐯[3][1] - 𝐯[2][1]) * (𝐱[2] - 𝐯[3][2])) / den
    𝐰[2] = ((𝐯[3][2] - 𝐯[1][2]) * (𝐱[1] - 𝐯[3][1]) +
             (𝐯[1][1] - 𝐯[3][1]) * (𝐱[2] - 𝐯[3][2])) / den
    𝐰[3] = 1 - 𝐰[1] - 𝐰[2]
    # end 计算重心坐标权重

    return 𝐰 .* 𝐁s # 用权重 𝐰 对 𝐁 插值
end

function discrete(x,y,z)
    R, ϕ, Z = cart2cyld(x,y,z)
    iϕ::Int = ϕ ÷ Δϕ
    dϕ = (ϕ-iϕ*Δϕ)/Δϕ
    𝐱 = [R-R0, Z]
    three = zeros(Int, 3)
    v_id = get_vertex_id_by_pos!(three, 𝐱, MD)
    
    if v_id == -1
        𝐁_prev = sum(interpolate𝐁(three, 𝐱, @view(B[:, :, iϕ])))
        𝐁_next = sum(interpolate𝐁(three, 𝐱, @view(B[:, :, cycle(iϕ+1, 18)])))
        
        if use_electric_field
            E_prev = sum(interpolate𝐁(three, 𝐱, @view(E[:, :, iϕ])))
            E_next = sum(interpolate𝐁(three, 𝐱, @view(E[:, :, cycle(iϕ+1, 18)])))
            return (parent((1-dϕ) * 𝐁_prev + dϕ * 𝐁_next), 
                   parent((1-dϕ) * E_prev + dϕ * E_next))
        else
            return (parent((1-dϕ) * 𝐁_prev + dϕ * 𝐁_next), 
                   zeros(3))  # 当不使用电场时返回零电场
        end
    end
    
    v_id == -2 && error("出界!")
    error("error!")
end

end # module
