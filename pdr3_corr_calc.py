import pandas as pd
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
import matplotlib.pyplot as plt
import numpy as np

input_path = "/mnt/hdd1/Projects/spring_gdrop_bao/"

data = pd.read_csv(input_path+"springfull_data.txt", names = ['ra', 'dec', 'mag'], sep = ' ')
#data = data[(data['ra'] >= 154) & (data['ra'] <= 225) & (data['dec'] >= -1.4) & (data['dec'] <= 1.0)]
data = data[(data['ra'] >= 160) & (data['ra'] <= 170) & (data['dec'] >= -1.4) & (data['dec'] <= 1.0)]
data_rand = pd.read_csv(input_path+"springfull_rand.txt", names = ['ra', 'dec'], sep = ' ')
#data_rand = data_rand[(data_rand['ra'] >= 154) & (data_rand['ra'] <= 225) & (data_rand['dec'] >= -1.4) & (data_rand['dec'] <= 1.0)]
data_rand = data_rand[(data_rand['ra'] >= 160) & (data_rand['ra'] <= 170) & (data_rand['dec'] >= -1.4) & (data_rand['dec'] <= 1.0)]
data_rand = data_rand.sample(frac = 0.1)

mag_bin1 = data[(data['mag'] > 24.0) & (data['mag'] <= 25.0)]
mag_bin2 = data[(data['mag'] > 25.0) & (data['mag'] <= 26.0)]

def rand_theta_pairs(rand_data, rbins):
    nthreads = 4
    autocorr = 1
    RR_theta_pairs = DDtheta_mocks(autocorr, nthreads, rbins, rand_data['ra'], rand_data['dec'], verbose = True)
    return RR_theta_pairs
    


def calc_w_theta(autocorr, data_list, rbins, RR_counts):
    # data_list should be provided as a list even if there is only element
    
    nthreads = 4
    
    if autocorr == 1:
        data1 = data_list[0]
        DDtheta_pairs = DDtheta_mocks(autocorr, nthreads, rbins, data1['ra'], data1['dec'], verbose = True)
        ang_corr = DDtheta_pairs['npairs']/(frac*frac*RR_counts) - 1
        ang_corr_err = np.sqrt(DDtheta_pairs['npairs'])/(frac*frac*RR_counts)

        return ang_corr, ang_corr_err

    if autocorr == 0:

        data1 = data_list[0]
        data2 = data_list[1]
        DDtheta_pairs = DDtheta_mocks(autocorr, nthreads, rbins, data1['ra'], data1['dec'], RA2 = data2['ra'], DEC2 = data2['dec'], verbose = True)
        ang_corr = DDtheta_pairs['npairs']/(frac*frac*RR_counts) - 1
        ang_corr_err = np.sqrt(DDtheta_pairs['npairs'])/(frac*frac*RR_counts)
        return ang_corr, ang_corr_err

frac = len(data)/len(data_rand)

rbins = np.linspace(0, 2, 21)
bin_mid = (rbins[1:] +rbins[0:-1])/2
RR_theta_pairs = rand_theta_pairs(data_rand, rbins)
corr_mag_bin1, corr_mag_bin1_err = calc_w_theta(1, [mag_bin1], rbins, RR_theta_pairs['npairs'])

#plt.scatter(data_rand['ra'], data_rand['dec'], s = 2)
#plt.savefig('pdr3_spring_crop_rand.png')

plt.errorbar(bin_mid, corr_mag_bin1, yerr = corr_mag_bin1_err, fmt = '.')
#plt.yscale('log')
plt.show()
