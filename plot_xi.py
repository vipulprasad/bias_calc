import numpy as np
import matplotlib.pyplot as plt 


def plot_ratio(mbin1, mbin2):

        file11 = '2pnt_cross_corr_M_{}_{}.txt'.format(mbin1, mbin1)
        file22 = '2pnt_cross_corr_M_{}_{}.txt'.format(mbin2, mbin2)
        file12 = '2pnt_cross_corr_M_{}_{}.txt'.format(mbin1, mbin2)

        print(file11+'\n', file22+'\n', file12+'\n')

        data11 = np.loadtxt(input_path+file11)
        distance11 = data11[:,0]
        corr_2pnt11 = data11[:,1]

        data22 = np.loadtxt(input_path+file22)
        distance22 = data22[:,0]
        corr_2pnt22 = data22[:,1]

        data12 = np.loadtxt(input_path+file12)
        distance12 = data12[:,0]
        corr_2pnt12 = data12[:,1]

        ratio = np.divide(np.sqrt(corr_2pnt11*corr_2pnt22), corr_2pnt12)

        plt.figure(figsize=(8, 6))
        plt.title('Mass bins: {}, {}'.format(mbin1, mbin2))
        plt.plot([0.1, 50], [1., 1.],'--',color='orange')
        plt.scatter(distance11, ratio, s = 16)#, yerr = corr_2pnt[:,2])
        #plt.errorbar(distance, corr_2pnt, yerr = corr_2pnt[:,2])
        plt.xscale('log')
        #plt.yscale('log')
        plt.ylim((0.8, 1.2))
        plt.xlim((0, 50))
        plt.xlabel('r')
        plt.ylabel(r'$\frac{\sqrt{\xi_{hh}(M1)\xi_{hh}(M2)}}{\xi_{hh}(M1,M2)}$')
        plt.savefig(output_path+'corr_ratio_{}_{}_.png'.format(mbin1, mbin2))
        plt.show()

input_path = '/home/vipul/vipul/halo_clutering/bias_calc/abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z0.300/'
output_path = '/home/vipul/vipul/halo_clutering/bias_calc/abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z0.300/'

bin1 = '1E+12_1.5E+12'
bin2 = '2E+12_3E+13' 

plot_ratio(bin1, bin2)