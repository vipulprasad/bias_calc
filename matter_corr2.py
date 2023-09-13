import numpy as np
import matplotlib.pyplot as plt 
from colossus.cosmology import cosmology

def matter_corr(distance, z):

        cosmo = cosmology.setCosmology('planck18')
        xi_mm = cosmo.correlationFunction(distance, z = z)
        return xi_mm

def bias_calc(file):
        file_names = open(bias_filenames_file, 'a')
        m1, m2 = file.split('_')[3], file.split('_')[4]
        data = np.loadtxt(input_path+file)
        #print(file)

        distance = data[:,0]
        two_pnt_func = data[:,1]
        matter_corr_func =  matter_corr(distance, redshift)
        bias = np.sqrt(np.divide(two_pnt_func, matter_corr_func))
        bias_err = np.divide(data[:,2], 2*np.sqrt(two_pnt_func*matter_corr_func))
        #bias_err = np.sqrt(data[:,2], matter_corr_func)
        bias_result = np.stack((distance, bias, bias_err), axis = 1)
        out_fname = "bias{}".format(file.split('M')[1])
        file_names.write(out_fname+'\n')
        file_names.close()
        print(out_fname)
        np.savetxt(output_path+out_fname, bias_result)
        #return bias, m1, m2
#------------------------------------------------------------------#

input_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/corr_files/"
output_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/bias_files/"

corr2pnt_filenames_file = "2pnt_file_names_7500.txt"
bias_filenames_file = 'bias_files_7500.txt'

corr2pnt_filenames = open(corr2pnt_filenames_file, 'r').readlines()
radii = np.logspace(-2, 3, 100) #Mpc/h comoving
redshift = 3.0

open(bias_filenames_file, 'w').close()

for filename in corr2pnt_filenames:
        bias_calc(filename.split('\n')[0])
