import numpy as np
import matplotlib.pyplot as plt 
from colossus.cosmology import cosmology


def matter_corr(distance, z):

        cosmo = cosmology.setCosmology('planck18')
        xi_mm = cosmo.correlationFunction(distance, z = z)
        return xi_mm


radii = np.logspace(-2, 3, 100) #Mpc/h comoving
#Assume that k and P come from somewhere, e.g. CAMB or CLASS
#xi_mm = xi.xi_mm_at_r(radii, kh, pk) 



def bias_calc(file):
        file_name = '/home/vipul/vipul/halo_clutering/bias_marked/bias_plots_0.01-30/corr_files/{}'.format(file)
        m1, m2 = file.split('_')[3], file.split('_')[4]
        data = np.loadtxt(file_name)
        #print(file)

        distance = data[:,0]
        two_pnt_func = data[:,1]
        matter_corr_func =  matter_corr(distance, 3.0)
        bias = np.sqrt(np.divide(two_pnt_func, matter_corr_func))
        bias_err = np.divide(data[:,2], 2*np.sqrt(two_pnt_func*matter_corr_func))
        #bias_err = np.sqrt(data[:,2], matter_corr_func)
        bias_result = np.stack((distance, bias, bias_err), axis = 1)
        out_fname = "bias_plots_0.01-30/bias_files/bias{}".format(file.split('\n')[0].split('M')[1])
        print(out_fname)
        np.savetxt(out_fname, bias_result)
        #return bias, m1, m2

#bias_calc('2pnt_corr_M_8.01E+10_8.65E+10.txt')
#bias_val, ml, mu = bias_calc('2pnt_corr_M_8.01E+10_8.65E+10.txt')
#bias_calc('2pnt_cross_corr_M_1.85E+11_2E+11_1.28E+11_1.38E+11.txt')

file_names = open('2pnt_file_names.txt', 'r').readlines()

for file in file_names:
        bias_calc(file.split('\n')[0])