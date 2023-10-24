import numpy as np
import matplotlib.pyplot as plt 
import argparse
import scienceplots
plt.style.use('science')

def plot(mbin):

        corr_file = '2pnt_cross_corr_M_{}_.txt'.format(mbin)
        data = np.loadtxt(input_path+file11)
        distance = data[:,0]
        corr_2pnt = data[:,1]
        corr_2pnt_error = data[:,2]

        #plt.scatter(distance, corr_2pnt, s = 8)#, yerr = corr_2pnt[:,2])
        plt.errorbar(distance, corr_2pnt, yerr = corr_2pnt_error, fmt = '.')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r')
        plt.ylabel(r'$\xi_{hh}(r)$')
        plt.show()

def plot_ratio(mbin1, mbin2):

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
        plt.vlines(exclsn, -1, 2, color = 'black', linestyles=':')
        #plt.scatter(distance11, ratio, s = 16)#, yerr = corr_2pnt[:,2])
        plt.errorbar(distance11, ratio, yerr = ratio_err, fmt = '.')
        plt.xscale('log')
        #plt.yscale('log')
        plt.ylim((0.8, 1.2))
        plt.xlim((0.1, 50))
        plt.xlabel('r', fontsize = 18)
        plt.ylabel(r'$\frac{\sqrt{\xi_{hh}(M1)\xi_{hh}(M2)}}{\xi_{hh}(M1,M2)}$', fontsize = 18)
        print(output_path+'corr_ratio_{}_{}_.png'.format(mbin1, mbin2))
        plt.savefig(output_path+'corr_ratio_{}_{}_.png'.format(mbin1, mbin2), format="png")
        plt.show()

parser = argparse.ArgumentParser()
parser.add_argument('z', type = float, help = 'redshift')
args = parser.parse_args()

redshift = args.z


if redshift in [0.3, 0.5]:
        input_path = 'abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z{}00/corr_files/'.format(redshift)
        output_path = 'abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z{}00/ratio_plots/'.format(redshift)
        sim_boxsize = 1100.0
        exclsn = 2.8827501554159545

if redshift in [3.0, 4.0]:
        input_path = 'abacus_summit/box7500/z{}/corr_files/'.format(redshift)
        output_path = 'abacus_summit/box7500/z{}/ratio_plots/'.format(redshift)
        sim_boxsize = 7500.0
        exclsn = 1.1127327321765554 #halo exclusion scale in MPc

massbin_filenames_file = 'massbin_z{}.txt'.format(redshift)

file_names = open(massbin_filenames_file, 'r').readlines()
file_names = [file.split('\n')[0] for file in file_names]

m_bins = []

print("\n")

for i, file_name in enumerate(file_names):
        print("({}) - {}".format(i, file_name))
        m_bins.append(file_name.split('_')[2])

mbin1 = m_bins[int(input("\nSelect first mass bin\n"))]
mbin2 = m_bins[int(input("\nSelect second mass bin\n"))]

print("First mass_bin - {}".format(mbin1))
print("Second mass_bin - {}".format(mbin2))

plot_ratio(mbin1, mbin2)
