module UserInputs
# struct Inputs
#     Δt::Float64
#     N::Int # number of particles
#     TotalSteps::Int # total steps
# end

use_electric_field = false  # 控制是否使用电场的开关
is_data_saving_on = true
is_merge_process_files_on = true

# 初始化类型
init_type = :single  # 或 :parabolic_torus

# 粒子数量（分布模式下使用）
N = 1000

# 动量分布参数
μ = 5.0        # 动量大小的平均值
σ = 2.0        # 动量大小的标准差
sinθ_min = 0.0 # 最小投掷角的正弦值
sinθ_max = 1.0 # 最大投掷角的正弦值
ϕ_min = 0.0    # 最小回旋角
ϕ_max = 2π     # 最大回旋角

# 其他现有配置...

include("user_inputs.jl")
# inputs=Inputs(Δt, N, TotalSteps)
end
