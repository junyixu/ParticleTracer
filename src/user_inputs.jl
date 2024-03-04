import PhysicalConstants.CODATA2018 as C
Δt = 1.0e-1
N = 1 # number of particles
TotalSteps = 100000 # total steps
B0 = 2.0 # Magnetic strength (T)
E0 = 2.0 # Electric strength (V/m)
R0 = 1.7 # Major radius of torus (m)
a  = 0.4 # Minor radius of torus (m)
x0 = [1.8, 0, 0] # initial position (m)
unit_p=C.m_e.val*C.c_0.val
p0 = [1, 5, 0]*unit_p # initial momentum (kg*m/s)
ptc_type="electron" # particle type
