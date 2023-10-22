import pandas as pd
import numpy as np
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse

def format_e(n):
    a = '%.2E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


#--------------------------auto correlation calculations----------------------------------#
nthreads = 4

def rand_npairs(rand_N, boxsize, bins):
        
        rand_X = np.random.uniform(0, boxsize, rand_N)
        rand_Y = np.random.uniform(0, boxsize, rand_N)
        rand_Z = np.random.uniform(0, boxsize, rand_N)

        autocorr=1
        RR_counts = DD(autocorr, nthreads, bins, rand_X, rand_Y, rand_Z, periodic = True, verbose = True, boxsize = boxsize)
        return RR_counts

def rand_npairs_2(N, boxsize, bins):

        volume = boxsize**3 
        pair_counts = []
        for i in range(len(bins) -1 ):
            pair_counts.append(N*(N-1)*(4/3)*np.pi*(bins[i+1]**3 - bins[i]**3)/volume)
        return np.array(pair_counts)

def calc_xi(file_names, autocorr, bins):

        RR_counts = rand_npairs_2(rand_N, sim_boxsize, bins)
        
        if autocorr == 1: 
                mbin1 = file_names[0].split('_')[2]
                print("mass bin edge = {} - {}".format(Ml, Mu))
                
                data = pd.read_csv(input_path+file_names[0], names = ["M", "X", "Y", "Z"])
                f = len(data)/rand_N

                print("mass bin size = ",len(data))
                data_X = data['X']
                data_Y = data['Y']
                data_Z = data['Z']

                auto_DD_counts = DD(autocorr, nthreads, bins, data_X, data_Y, data_Z, periodic = True, boxsize = sim_boxsize, verbose = True)
                del data, data_X, data_Y, data_Z

                auto_corr_2pnt = (auto_DD_counts['npairs']/((f**2)*RR_counts))-1
                #auto_corr_2pnt = (auto_DD_counts['npairs']/((f**2)*RR_counts['npairs']))-1
                auto_corr_2pnt_err = np.sqrt(auto_DD_counts['npairs'])/((f**2)*RR_counts)
                #auto_corr_2pnt_err = np.sqrt(auto_DD_counts['npairs'])/((f**2)*RR_counts['npairs'])
                auto_corr_2pnt = np.stack((r_bin_mid, auto_corr_2pnt, auto_corr_2pnt_err), axis = 1)

                out_fname = "2pnt_corr_M_{}_.txt".format(format_e(mbin1))
                np.savetxt(output_path+out_fname, auto_corr_2pnt)

                return auto_corr_2pnt

        if autocorr == 0:
                mbin1, mbin2 = file_names[0].split('_')[2], file_names[1].split('_')[2]
                
                data1 = pd.read_csv(input_path+file_names[0], names = ["M", "X", "Y", "Z"])
                data2 = pd.read_csv(input_path+file_names[1], names = ["M", "X", "Y", "Z"])
                f = len(data1)*len(data2)/rand_N**2

                data_X1 = data1['X']
                data_Y1 = data1['Y']
                data_Z1 = data1['Z']
                data_X2 = data2['X']
                data_Y2 = data2['Y']
                data_Z2 = data2['Z']

                auto_DD_counts = DD(autocorr, nthreads, bins, data_X1, data_Y1, data_Z1, X2 = data_X2, Y2 = data_Y2, Z2 = data_Z2, periodic = True, boxsize = sim_boxsize, verbose = True)
                del data1, data2, data_X1, data_Y1, data_Z1, data_X2, data_Y2, data_Z2

                #cross_corr_2pnt = (auto_DD_counts['npairs']/((f)*RR_counts['npairs']))-1
                #cross_corr_2pnt_err = np.sqrt(auto_DD_counts['npairs'])/(f*RR_counts['npairs'])
                cross_corr_2pnt = (auto_DD_counts['npairs']/((f)*RR_counts))-1
                cross_corr_2pnt_err = np.sqrt(auto_DD_counts['npairs'])/(f*RR_counts)
                cross_corr_2pnt = np.stack((r_bin_mid, cross_corr_2pnt, cross_corr_2pnt_err), axis = 1)

                out_fname = "2pnt_cross_corr_M_{}_{}_.txt".format(mbin1, mbin2)
                print(out_fname)
                np.savetxt(output_path+out_fname, cross_corr_2pnt)
                return cross_corr_2pnt
                

#------Matter correlation----------------- 


def plot(corr_2pnt):
        plt.scatter(bin_mid, corr_2pnt[:,1], s = 8)#, yerr = corr_2pnt[:,2])
        plt.errorbar(bin_mid, corr_2pnt[:,1], yerr = corr_2pnt[:,2])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r')
        plt.ylabel(r'$\xi_{hh}(r)$')
        plt.show()

def do_calc_xi(file_names):

        for file1 in file_names:
                for file2 in file_names:
                        calc_xi([file1, file2], autocorr =0, bins = r_bins)

def do_calc_xi_pool(file_names):

        #open(corr2pnt_filenames_file, "w").close()
        #file_names = open(massbin_filenames_file, 'r').readlines()
        #file_names = [file.split('\n')[0] for file in file_names]

        cross_corr_input_arg = []

        for file1 in file_names:
            for file2 in file_names:
                cross_corr_input_arg.append(([file1, file2], 0))

        with Pool(3) as P:
                P.starmap(calc_xi, cross_corr_input_arg)        


if __name__ == "__main__":

        parser = argparse.ArgumentParser()
        parser.add_argument('z', type = float, help = 'redshift')
        args = parser.parse_args()

        redshift = args.z

        if redshift in [0.3, 0.5]:
                input_path = 'abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z{}00/'.format(redshift)
                output_path = 'abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z{}00/corr_files/'.format(redshift)
                sim_boxsize = 1100.0

        if redshift in [3.0, 4.0]:
                input_path = 'abacus_summit/box7500/z{}/halo_cat/'.format(redshift)
                output_path = 'abacus_summit/box7500/z{}/corr_files/'.format(redshift)
                sim_boxsize = 7500.0

        massbin_filenames_file = 'massbin_z{}.txt'.format(redshift)

        file_names = open(massbin_filenames_file, 'r').readlines()
        file_names = [file.split('\n')[0] for file in file_names]

        rand_N = int(1e6)
        nbins = 20
        r_bins = np.logspace(-1, 1.5, nbins +1)
        r_bin_mid = (r_bins[0:-1] + r_bins[1:])/2

        do_calc_xi(file_names)


