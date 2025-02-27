module UserInputs
# struct Inputs
#     Δt::Float64
#     N::Int # number of particles
#     TotalSteps::Int # total steps
# end

use_electric_field = false  # 控制是否使用电场的开关
include("user_inputs.jl")
# inputs=Inputs(Δt, N, TotalSteps)

end
