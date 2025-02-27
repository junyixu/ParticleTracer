# simulation parameters
Δt = 1.0
N = 1 # number of particles
# TotalSteps = ceil(Int, 40000*2π) # total steps
TotalSteps = 10
SavePerNSteps = 1 # Save 1000 steps

# tokamak parameters
B0 = 2.0 # Magnetic strength (T)
E0 = 2.0 # Electric strength (V/m)
R0 = 1.7 # Major radius of torus (m)
a  = 0.4 # Minor radius of torus (m)

u = Unit(B0)
B0 /= u.B # Magnetic strength (T)
E0 /= u.E # Electric strength (V/m)
R0 /= u.x # Major radius of torus (m)
a  /= u.x # Minor radius of torus (m)

# initial conditions
init_type = :single
x0 = Float64[1.8, 0, 0] / u.x # initial position (m)
p0 = Float64[5.0, 1, 0] # initial momentum (kg*m/s)

ptc_type=:electron # particle type

# pusher=:boris # pusher type
# pusher=:RVPA_Cay3D
pusher=:relativistic_boris_step
