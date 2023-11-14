from math import pi, pow


def calc_exclusion(halo_mass):

        Mpc = 3.0857e22                     # 1 Mpc in metres.
        sm = 1.98e30                           # Solar mass in Kg.
        G = 6.67408e-11#*sm/(Mpc**3)         # Gravitational constant.
        h = 0.678                   # Hubble constant in units of 100 Km/s/Mpc.
        H_o = (100*h)/(Mpc*(10**(-3)))      # Hubble constant.
        rho_crit = (3*(H_o)**2)/(8*pi*G)*(Mpc**3/sm) # critical density of the universe.
        delta = 200

        scale = pow((3*halo_mass)/(4*pi*delta*rho_crit), 1/3)

        print("Halo exclusion scale = {} MPc".format(scale))

        return scale
