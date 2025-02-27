# simulation parameters
Δt = 0.1
N = 1 # number of particles
# TotalSteps = ceil(Int, 40000*2π) # total steps
TotalSteps = 100000
SavePerNSteps = 1 # Save 1000 steps

# tokamak parameters
B0 = 5.18 # Magnetic strength (T)
E0 = 2.0 # Electric strength (V/m)
R0 = 6.2 # Major radius of torus (m)
a  = 2.0 # Minor radius of torus (m)

u = Unit(B0, "alpha")
B0 /= u.B # Magnetic strength (T)
E0 /= u.E # Electric strength (V/m)
R0 /= u.x # Major radius of torus (m)
a  /= u.x # Minor radius of torus (m)

# initial conditions
x0 = Float64[1.57941 , 5.24294 , -1.29892] # initial position (m)
p0 = Float64[ -0.0284313 , -0.0247643 ,   0.0213838 ] # initial momentum (kg*m/s)

use_electric_field = true  # 控制是否使用电场的开关
ptc_type=:alpha # particle type

# pusher=:boris # pusher type
pusher=:RVPA_Cay3D
# pusher=:relativistic_boris_step
field=:discrete
