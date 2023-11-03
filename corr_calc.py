import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Corrfunc.theory.DD import DD
#from Corrfunc.utils import convert_3d_counts_to_cf
from multiprocessing import Pool
import argparse
from math import pi, pow
nthreads = 32

def format_e(n):
    a = '%.2E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


#--------------------------auto correlation calculations----------------------------------#

def halo_vir_rad(halo_mass): # function to find halo virial radius

        Mpc = 3.0857e22                     # 1 Mpc in metres.
        sm = 1.98e30                           # Solar mass in Kg.
        G = 6.67408e-11#*sm/(Mpc**3)         # Gravitational constant.
        h = 1#0.678                   # Hubble constant in units of 100 Km/s/Mpc.
        H_o = (100*h)/(Mpc*(10**(-3)))      # Hubble constant.
        rho_crit = (3*(H_o)**2)/(8*pi*G)*(Mpc**3/sm) # critical density of the universe.
        delta = 200

        scale = pow((3*halo_mass)/(4*pi*delta*rho_crit), 1/3)

        print("Halo exclusion scale = {} MPc".format(scale))

        return scale

def rand_npairs(rand_N, boxsize, bins):
        
        rand_X = np.random.uniform(0, boxsize, rand_N)
        rand_Y = np.random.uniform(0, boxsize, rand_N)
        rand_Z = np.random.uniform(0, boxsize, rand_N)

        autocorr=1
        RR_counts = DD(autocorr, nthreads, bins, rand_X, rand_Y, rand_Z, periodic = True, verbose = True, boxsize = boxsize)
        return RR_counts

def rand_npairs_2(N, boxsize, bins): # theoretical calculation or random-random pair counts

        volume = boxsize**3 
        pair_counts = []
        for i in range(len(bins) -1 ):
            pair_counts.append(N*(N-1)*(4/3)*np.pi*(bins[i+1]**3 - bins[i]**3)/volume)
        return np.array(pair_counts)

def calc_xi(file_name1, file_name2):

        mbin1, mbin2 = file_name1.split('_')[2], file_name2.split('_')[2]

        rand_N = int(1e6)
        nbins = 20
        r_bins = np.logspace(-1, 1.5, nbins +1)
        r_bin_mid = (r_bins[0:-1] + r_bins[1:])/2

        RR_counts = rand_npairs_2(rand_N, sim_boxsize, r_bins)
        
        data1 = pd.read_csv(input_path+file_names[0], names = ["M", "X", "Y", "Z"])
        data2 = pd.read_csv(input_path+file_names[1], names = ["M", "X", "Y", "Z"])
        f = len(data1)*len(data2)/rand_N**2

        data_X1 = data1['X']
        data_Y1 = data1['Y']
        data_Z1 = data1['Z']
        data_X2 = data2['X']
        data_Y2 = data2['Y']
        data_Z2 = data2['Z']

        autocorr = 0 # For cross correlation calculation. To calculate autocorrelation make file1 and file2 same.
        auto_DD_counts = DD(autocorr, nthreads, r_bins, data_X1, data_Y1, data_Z1, X2 = data_X2, Y2 = data_Y2, Z2 = data_Z2, periodic = True, boxsize = sim_boxsize, verbose = True)
        del data1, data2, data_X1, data_Y1, data_Z1, data_X2, data_Y2, data_Z2

        #cross_corr_2pnt = (auto_DD_counts['npairs']/((f)*RR_counts['npairs']))-1
        #cross_corr_2pnt_err = np.sqrt(auto_DD_counts['npairs'])/(f*RR_counts['npairs'])
        cross_corr_2pnt = (auto_DD_counts['npairs']/((f)*RR_counts))-1
        cross_corr_2pnt_err = np.sqrt(auto_DD_counts['npairs'])/(f*RR_counts)
        cross_corr_2pnt = np.stack((r_bin_mid, cross_corr_2pnt, cross_corr_2pnt_err), axis = 1)

        out_fname = "2pnt_cross_corr_M_{}_{}_.txt".format(mbin1, mbin2)
        print(out_fname)
        np.savetxt(output_path_corr+out_fname, cross_corr_2pnt)
        return None #cross_corr_2pnt

def calc_ratio(mbin1, mbin2):
        
        # halo exclusion scale is taken as the sum of the virial radius of largest halo of the two bins.
        exclusion = halo_vir_rad(float(mbin1.split('-')[1])) + halo_vir_rad(float(mbin2.split('-')[1]))


        file11 = '2pnt_cross_corr_M_{}_{}_.txt'.format(mbin1, mbin1)
        file22 = '2pnt_cross_corr_M_{}_{}_.txt'.format(mbin2, mbin2)
        file12 = '2pnt_cross_corr_M_{}_{}_.txt'.format(mbin1, mbin2)

        print(file11+'\n', file22+'\n', file12+'\n')

        data11 = np.loadtxt(input_path+file11)
        distance11 = data11[:,0]
        corr_2pnt11 = data11[:,1]
        corr_2pnt_err11 = data11[:,2]

        data22 = np.loadtxt(input_path+file22)
        distance22 = data22[:,0]
        corr_2pnt22 = data22[:,1]
        corr_2pnt_err22 = data22[:,2]

        data12 = np.loadtxt(input_path+file12)
        distance12 = data12[:,0]
        corr_2pnt12 = data12[:,1]
        corr_2pnt_err12 = data12[:,2]

        ratio = np.divide(np.sqrt(corr_2pnt11*corr_2pnt22), corr_2pnt12)
        ratio_err = ratio*((0.5*corr_2pnt_err11/corr_2pnt11) + (0.5*corr_2pnt_err22/corr_2pnt22) + (corr_2pnt_err12/corr_2pnt12))


        plt.figure(figsize=(10, 6))
        plt.title('Mass bins: {}, {}, z = {}'.format(mbin1, mbin2, redshift), fontsize = 20)
        plt.plot([0.1, 50], [1., 1.],'--',color='orange')
        plt.vlines(exclusion, -2,max(ratio)+1.0 , color = 'black', linestyles=':', label = 'halo exclusion scale= {exclusion} MPc')
        #plt.scatter(distance11, ratio, s = 16)#, yerr = corr_2pnt[:,2])
        plt.errorbar(distance11, ratio, yerr = ratio_err, fmt = '.')
        plt.xscale('log')
        #plt.yscale('log')
        plt.ylim((0.6, max(ratio)+0.5))
        plt.xlim((exclusion-0.5, 50))
        plt.xlabel('r', fontsize = 18)
        plt.ylabel(r'$\frac{\sqrt{\xi_{hh}(M1)\xi_{hh}(M2)}}{\xi_{hh}(M1,M2)}$', fontsize = 18)
        print(output_path_ratio+'corr_ratio_{}_{}_.png'.format(mbin1, mbin2))
        plt.savefig(output_path_ratio+'corr_ratio_{}_{}_.png'.format(mbin1, mbin2), format="png")
        plt.close()
                

def do_calc_xi(file_names):

        for i in range(len(file_names)):
                for j in range(len(file_names[i:])):
                        #calc_xi([file_names[i], file_names[i:][j]], bins = r_bins)
                        calc_xi(file_names[i], file_names[i:][j])

def do_calc_xi_pool(file_names):

        #open(corr2pnt_filenames_file, "w").close()
        #file_names = open(massbin_filenames_file, 'r').readlines()
        #file_names = [file.split('\n')[0] for file in file_names]

        cross_corr_input_arg = []

        for i in range(len(file_names)):
                for j in range(len(file_names[i:])):
                        cross_corr_input_arg.append((file_names[i], file_names[i:][j], 0))

        with Pool(3) as P:
                P.starmap(calc_xi, cross_corr_input_arg)      

def do_plot_all(m_bins):

        for i in range(len(m_bins)):
                for j in range(len(m_bins[i:])):
                        #calc_xi([file_names[i], file_names[i:][j]], bins = r_bins)
                        calc_ratio(m_bins[i], m_bins[i:][j])

if __name__ == "__main__":

        parser = argparse.ArgumentParser()
        parser.add_argument('z', type = float, help = 'redshift')
        parser.add_argument('box_type', type = str, help = 'base, huge, small, hugebase' )
        args = parser.parse_args()

        box_size_dict = {"small": 500.0, "base": 2000.0, "hugebase": 2000.0, "huge": 7500}
        part_mass_dict = {"small": 2.109081520453063e+09, "base": 2.109081520453063e+09, "hugebase": 5.694520105223270e+10, "huge": 5.694520105223272e+10} #hMsun
        #part_mass_dict = {"small": 3.131059264330557e+09, "base": 3.131059264330557e+09, "hugebase": 8.453860013692505e+10, "huge": 8.453860013692505e+10} #Msun

        h = 1#0.678
        redshift = args.z
        box_type = args.box_type
        sim_boxsize = box_size_dict[box_type]
        part_mass = part_mass_dict[box_type]

        input_path = 'z{}/{}/halo_cat/'.format(redshift, box_type)
        output_path_corr = 'z{}/{}/corr_files/'.format(redshift, box_type)
        output_path_ratio = 'z{}/{}/ratio_plots/'.format(redshift, box_type)

        massbin_filenames_file = 'massbin_z{}_{}.txt'.format(redshift, box_type)

        file_names = open(massbin_filenames_file, 'r').readlines()
        file_names = [file.split('\n')[0] for file in file_names]

        m_bins = []

        for i, file_name in enumerate(file_names):
                #print("({}) - {}".format(i, file_name))
                m_bins.append(file_name.split('_')[2])

        # rand_N = int(1e6)
        # nbins = 20
        # r_bins = np.logspace(np.log10(exclsn), 1.5, nbins +1)
        # r_bin_mid = (r_bins[0:-1] + r_bins[1:])/2

        do_calc_xi(file_names)


