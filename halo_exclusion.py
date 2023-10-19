from math import pi, pow


delta = 200

halo_mass = 147886687132648.38 # mass of smallest halo solar mass

omega_m = 0.301685                  # Matter(baryons+dark matter) density parameter.
Mpc = 3.0857e22                     # 1 Mpc in metres.
sm = 2e30                           # Solar mass in Kg.
G = 6.67408e-11*sm/(Mpc**3)         # Gravitational constant.
h = 0.674                   # Hubble constant in units of 100 Km/s/Mpc.
H_o = (100*h)/(Mpc*(10**(-3)))      # Hubble constant.
rho_crit = (3*(H_o)**2)/(8*pi*G) # critical density of the universe.
rho_bar = omega_m*rho_crit          # comoving background density in kg/m^3.


scale = pow((3*halo_mass*sm)/(4*pi*delta*rho_crit), 1/3)
print(scale)