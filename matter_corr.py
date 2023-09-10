import numpy as np
from math import sin, cos, exp, pi, sqrt, log
import matplotlib.pyplot as plt
import argparse

Mpc = 3.0857e22                     # 1 Mpc in metres.
sm = 2e30                           # Solar mass in Kg.
G = 6.67408e-11*sm/(Mpc**3)         # Gravitational constant.
h = 0.674                           # Hubble constant in units of 100 Km/hr/Mpc.
H_o = (100*h)/(Mpc*(10**(-3)))      # Hubble constant.
omega_b = 0.0224/(h**2)             # Baryons density parameter.
omega_m = 0.315                     # Matter(baryons+dark matter) density parameter.
omega_lambda = 0.6847               # Dark energy density parameter.
deltac = 1.686                      # Critical overdensity for spherical collapse at z = 0.
sigma_8 = 0.811                     # sigma_8 used for normalization.
ns = 0.965                          # spectral indece of initial power spectrum.
rho_crit = (3*(H_o)**2)/(8*np.pi*G) # critical density of the universe.
rho_bar = omega_m*rho_crit          # comoving background density in kg.


def Simpson(f, a, b, m): #if a = 0 this won't work
	x = a				 #also for small intervals Simpson has more relative 
	h  = a*1.01			 #error than Trapezoidal below
	s = 0
	
	while x<b:
		s += (h/3)*(f(x,m)+4*f(x+h,m)+f(x+2*h,m))
		x += 2*h
		h*=1.01
		
	return s
	
def Trapezoidal(f, a, b, m): 
	#if a = 0 this won't work
	x = a
	h  = a*1.01
	s = 0
	while x<b:
		s += (h/2)*(f(x,m)+f(x+h,m))
		x += h
		h*=1.01
		
	return s


def Transfer(k): # Transfer function

	gamma = h*h*omega_m*exp(-omega_b-(sqrt(2*h)*(omega_b/omega_m)))
	x = (k/gamma)
	A = log(1+2.34*x)/(2.34*x)
	B = (1+(3.89*x)+(16.1*x)**2+(5.46*x)**3+(6.71*x)**4)**(-0.25)
	return A*B

def D(z): 
	# Growth factor normalized to unity at z = 0.
	
	A = omega_m+0.4545*omega_lambda
	B = omega_m*(1+z)**3 + 0.4545*omega_lambda
	return (A/B)**(1/3)

def deltac(z): 
	# Critical overdensity for spherical collapse at different epoch.
	return 1.686/D(z)

def Power_spec(k):
	# Power spectrum
	
	P_i = k**ns # Initial power spectrum
	T = Transfer(k)**2
	return (norm8**2)*P_i*T 

def R(m):
	# comoving radius of a sphere of mass m
	return ((3*m)/(4*PI*rho_bar))**(1/3)

def W_kR(k,R):
	# window function in fourier space, used to smoothen the density field
	# This particular one is the Fourier transform of tophat filter function.
	
	x = k*R
	A = sin(x)
	B = cos(x)
	return (3/(x**3))*(A-x*B)

def sigma_integrand(k,R):
	A = (W_kR(k,R))**2.0
	B = (k**2)*Power_spec(k)
	return A*B/(2*pi**2)

def Sigma(R):
	integral = Simpson(sigma_integrand, 1e-7, 1e7, R)
	return sqrt(integral)

def corr_2pnt_integrand(k, r):
       return k*Power_spec(k)*np.sin(k*r)/r 

def corr_2pnt(r, z):
       return (D(z)**2)*Simpson(corr_2pnt_integrand, 1e-7, 1e7, r)/(2*np.pi**2)

def bias_calc(file):
        file_name = '/home/vipul/vipul/halo_clutering/bias_marked/{}'.format(file)
        m1, m2 = file.split('_')[3], file.split('_')[4]
        data = np.loadtxt(file_name)
        print(data[0:5])

        distance = data[:,0]
        two_pnt_func = data[:,1]
        matter_corr_func = [corr_2pnt(i, z) for i in distance]
        print(matter_corr_func)
        bias = np.sqrt(np.divide(two_pnt_func, matter_corr_func))


        plt.scatter(distance, bias, s = 8)
        plt.title('mass bin [{} - {}]'.format(m1, m2))    
        #plt.xscale('log')
        #plt.yscale('log')
        plt.xlabel(r'$r$')
        plt.ylabel(r'$b(r)$')
        plt.show()
        return bias, m1, m2

if __name__ == "__main__":

        z = 3
        norm8 = 1
        R8 = 8.0/h
        sigma8 = Sigma(R8)
        norm8= sigma_8/sigma8
        print(R8,Sigma(R8), norm8) 

        #bias_calc('2pnt_cross_corr_M_1.19E+11_1.28E+11_1.19E+11_1.28E+11.txt')
        #bias_calc('2pnt_corr_M_8.01E+10_8.65E+10.txt')
        matter_corr_my = [corr_2pnt(i, z) for i in distance]

