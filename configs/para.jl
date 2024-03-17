# simulation parameters
Δt = 1.0e-1
N = 8 # number of particles
# TotalSteps = 4000000 # total steps
# SavePerNSteps = 10000 # Save 1000 steps
TotalSteps = 100 # total steps
SavePerNSteps = 1 # Save 1000 steps

# tokamak parameters
B0 = 2.0 # Magnetic strength (T)
u = Unit(B0)
B0 = B0 * u.B
E0 = 2.0 * u.E # Electric strength (V/m)
R0 = 1.7 * u.x # Major radius of torus (m)
a  = 0.4 * u.x # Minor radius of torus (m)

# initial conditions

ptc_type=:electron # particle type

# pusher=:boris # pusher type
pusher=:RVPA_Cay3D

μ = 5.0
σ = 2.0
sinθ_min = 0.0
sinθ_max = 1.0
ϕ_min=0.0
ϕ_max=2π
