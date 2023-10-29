from math import pi, pow


delta = 200

halo_mass = 147886687132648.38 # mass of smallest halo solar mass 1.1127327321765554 MPc
#halo_mass = 2571446556156185.5 # z = 0.3 2.8827501554159545 MPc

omega_m = 0.301685                  # Matter(baryons+dark matter) density parameter.
Mpc = 3.0857e22                     # 1 Mpc in metres.
sm = 1.98e30                           # Solar mass in Kg.
G = 6.67408e-11#*sm/(Mpc**3)         # Gravitational constant.
h = 0.678                   # Hubble constant in units of 100 Km/s/Mpc.
H_o = (100*h)/(Mpc*(10**(-3)))      # Hubble constant.
rho_crit = (3*(H_o)**2)/(8*pi*G)*(Mpc**3/sm) # critical density of the universe.
rho_bar = omega_m*rho_crit          # comoving background density in kg/m^3.


def calc_exclusion(halo_mass):
        scale = pow((3*halo_mass)/(4*pi*delta*rho_crit), 1/3)

        if __name__ == "__main__":
                print("Halo exclusion scale = {} MPc".format(scale))

        return scale
