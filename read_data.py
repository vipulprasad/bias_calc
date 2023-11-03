import numpy as np
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import matplotlib.pyplot as plt
import argparse
plt.style.use('science')
import os

def format_e(n):
    a = '%.2E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def readdata(file_name, clms):

        print("Reading the data")

        cat = CompaSOHaloCatalog(file_name,cleaned=False)
        data = cat.halos[clms]
        data['N'] = data['N'].astype('float')*part_mass
        data['SO_central_particle'] += (box_size/2) # correcting the box coordinates origin
        data['X'] =data['SO_central_particle'][:,0]
        data['Y'] =data['SO_central_particle'][:,1]
        data['Z'] =data['SO_central_particle'][:,2]
        del cat
        return data[['N', 'X', 'Y', 'Z']].to_pandas()

def data_summary(data):

        bins = np.logspace(np.log10(min(data['N'])), np.log10(max(data['N'])), 100)
        bin_mid = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        plt.figure(figsize = (8,8))
        num, _ = np.histogram(data['N'], bins = bins)
        plt.scatter(bin_mid, 34*num)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("M")
        plt.ylabel("Number of halos")
        plt.show()

        print('tot particles = ' , len(data))
        print('min mass =', format_e(min(data['N'])))
        print('max mass =', format_e(max(data['N'])))

def create_mass_slice(data, ml_cut, mu_cut):

        outfile_file = open(massbin_filenames_file, 'a')

        if rand_sample != 0:
                rand_index = np.random.randint(0, len(data), rand_sample)
                data = data[rand_index]
                out_fname = "halos_mcut_{}-{}_{}_sample.txt".format(format_e(ml_cut), format_e(mu_cut), format_e(rand_sample))
                outfile_file.write(out_fname+'\n')
                print(out_fname)
        else:
                out_fname = "halos_mcut_{}-{}_.txt".format(format_e(ml_cut), format_e(mu_cut))
                outfile_file.write(out_fname+'\n')
                print(out_fname)

        outfile_file.close()

        data = data[(data['N'] > ml_cut) & (data['N'] <= mu_cut)]

        data.to_csv(output_path+out_fname, sep = ',', mode = 'a', index = False, header = False)

def load_files(mass_list):
        # To do -> Create seperate mass list
        halo_files = os.listdir(input_path)
        
        for file in halo_files:
                if file[0:4] == 'halo':
                        open(massbin_filenames_file, 'w').close()
                        print(file)
                        halo_cat = readdata(input_path+file, ['N', 'SO_central_particle'])
                        
                        for i in range(len(mass_list)-1):
                                create_mass_slice(halo_cat, mass_list[i], mass_list[i+1])
                        del halo_cat

#---------------------------------------------------------------------------#

parser = argparse.ArgumentParser()
parser.add_argument('z', type = float, help = 'redshift')
parser.add_argument('box_type', type = str, help = 'base, huge, small, huge_base' )
args = parser.parse_args()

box_size_dict = {"small": 500.0, "base": 2000.0, "hugebase": 2000.0, "huge": 7500}
part_mass_dict = {"small": 2.109081520453063e+09, "base": 2.109081520453063e+09, "hugebase": 5.694520105223270e+10, "huge": 5.694520105223272e+10} #hMsun
#part_mass_dict = {"small": 3.131059264330557e+09, "base": 3.131059264330557e+09, "hugebase": 8.453860013692505e+10, "huge": 8.453860013692505e+10} #Msun

h = 1#0.674
redshift = args.z
box_type = args.box_type
sim_boxsize = box_size_dict[box_type]
part_mass = part_mass_dict[box_type]

input_path = '/mnt/data/DATA/Simulations/AbacusSummit_Public_Data_Access/AbacusSummit_{}_c000_ph000/halos/z{}00/halo_info/'.format(box_type, redshift)
output_path = 'z{}/{}/halo_cat/'.format(redshift, box_type)

massbin_filenames_file = 'massbin_z{}_{}.txt'.format(redshift, box_type)

mass_list = [2e12, 2.2e12, 2.4e12, 2.6e12, 2.8e12, 3.2e12, 4e12, 6e12, 10e12]
#mass_list = [(3e12)*(5**i) for i in range(3)]+[max_halo_mass]
load_files(mass_list)


