import pandas as pd
import numpy as np
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf
import matplotlib.pyplot as plt
from multiprocessing import Pool

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

def calc_xi(file_names, autocorr):

        RR_counts = rand_npairs_2(rand_N, sim_boxsize, bins)

        outfile_file = open(corr2pnt_filenames_file, 'a')
        
        if autocorr == 1: 
                Ml, Mu = float(file_names[0].split('_')[2]), float(file_names[0].split('_')[3])
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
                auto_corr_2pnt = np.stack((bin_mid, auto_corr_2pnt, auto_corr_2pnt_err), axis = 1)
                out_fname = "2pnt_corr_M_{}_{}.txt".format(format_e(Ml), format_e(Mu))
                outfile_file.write(out_fname+'\n')
                np.savetxt(output_path+out_fname, auto_corr_2pnt)
                return auto_corr_2pnt

        if autocorr == 0:
                Ml1, Mu1, Ml2, Mu2 = float(file_names[0].split('_')[2]), float(file_names[0].split('_')[3]), float(file_names[1].split('_')[2]), float(file_names[1].split('_')[3])
                
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
                cross_corr_2pnt = np.stack((bin_mid, cross_corr_2pnt, cross_corr_2pnt_err), axis = 1)
                out_fname = "2pnt_cross_corr_M_{}_{}_{}_{}.txt".format(format_e(Ml1), format_e(Mu1), format_e(Ml2), format_e(Mu2))
                outfile_file.write(out_fname+'\n')
                print(out_fname)
                np.savetxt(output_path+out_fname, cross_corr_2pnt)
                return cross_corr_2pnt
        outfile_file.close()
                

#------Matter correlation----------------- 


def plot(corr_2pnt):
        plt.scatter(bin_mid, corr_2pnt[:,1], s = 8)#, yerr = corr_2pnt[:,2])
        plt.errorbar(bin_mid, corr_2pnt[:,1], yerr = corr_2pnt[:,2])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r')
        plt.ylabel(r'$\xi_{hh}(r)$')
        plt.show()

def abacus_cosmo_calc():

        halo_files = ['halos_mcut_1E+12_1.5E+12_.txt', 'halos_mcut_2E+12_3E+13_.txt', 'halos_mcut_3E+13_2E+15_.txt']

        for file1 in halo_files:
                for file2 in halo_files:
                        calc_xi([file1, file2], autocorr =0)

def abaus_summit_calc():

        #open(corr2pnt_filenames_file, "w").close()
        #file_names = open(massbin_filenames_file, 'r').readlines()
        #file_names = [file.split('\n')[0] for file in file_names]

        file_names = []
        cross_corr_input_arg = []

        for file1 in file_names:
            for file2 in file_names:
                cross_corr_input_arg.append(([file1, file2], 0))

        with Pool(3) as P:
                P.starmap(calc_xi, cross_corr_input_arg)        


if __name__ == "__main__":

        massbin_filenames_file = 'file_names_cosmos_z3.txt'
        corr2pnt_filenames_file = "2pnt_file_names_cosmo_z3.txt"

        nbins = 20
        #bins = np.linspace(0.1, 20.0, nbins + 1) # note that +1 to nbins
        bins = np.logspace(-1, 1.5, nbins +1)
        bin_mid = (bins[0:-1] + bins[1:])/2
        sim_boxsize = 1100
        rand_N = int(1e6)

        #input_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/halo_cat/"
        #output_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/corr_files/"

        input_path = '/home/vipul/vipul/halo_clutering/bias_calc/abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z0.300/'
        output_path = '/home/vipul/vipul/halo_clutering/bias_calc/abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z0.300/'
        abacus_cosmo_calc()

'''
        open(corr2pnt_filenames_file, "w").close()

        file_names = open(massbin_filenames_file, 'r').readlines()
        file_names = [file.split('\n')[0] for file in file_names]
        cross_corr_input_arg = []

        for file1 in file_names:
            for file2 in file_names:
                cross_corr_input_arg.append(([file1, file2], 0))

        with Pool(3) as P:
                P.starmap(calc_xi, cross_corr_input_arg)
'''
        
        