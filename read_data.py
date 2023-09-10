import numpy as np
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')
import os

h = 1#0.674

#part_mass = 3.131059264330557e+09*h
#part_mass = 2.109081520453063e+09 #hmsun 500
part_mass = 5.694520105223272e+10 #hmsun 7500
box_size = 500 # hmpc
#ml_cut = 2e11
#mu_cut = 2.2e11
rand_sample = 0 #int(1e5)
output_path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/halo_cat/"


def format_e(n):
    a = '%.2E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

#---------------------Read data----------------------------------------#

#file_path = "/home/vipul/vipul/halo_clutering/correlation_calc/z3.000/halo_info/halo_info_000.asdf"
file_path = "/home/vipul/vipul/halo_clutering/AbacusSummit_huge_c000_ph201/halos/z3.000/halo_info/"

def readdata(file_name, clms):

        print("Reading the data")

        cat = CompaSOHaloCatalog(file_name,cleaned=False)
        data = cat.halos[clms]
        data['N'] = data['N'].astype('float')*part_mass
        data['SO_central_particle'] += 3750 # correcting the box coordinates origin
        data['X'] =data['SO_central_particle'][:,0]
        data['Y'] =data['SO_central_particle'][:,1]
        data['Z'] =data['SO_central_particle'][:,2]
        del cat
        return data[['N', 'X', 'Y', 'Z']].to_pandas()

def data_summary(data):

        bins = np.logspace(np.log10(min(data['N'])), np.log10(max(data['N'])), 100)
        plt.figure(figsize = (8,8))
        plt.hist(data['N'], bins = bins)
        #plt.xscale("log")
        #plt.yscale("log")
        plt.xlabel("M")
        plt.ylabel("Number of halos")
        plt.show()

        print('tot particles = ' , len(data))
        print('min mass =', format_e(min(data['N'])))
        print('max mass =', format_e(max(data['N'])))

def create_mass_slice(data, ml_cut, mu_cut):

        open('file_names_7500.txt', 'w').close()
        outfile_file = open('file_names_7500.txt', 'a')

        if rand_sample != 0:
                rand_index = np.random.randint(0, len(data), rand_sample)
                data = data[rand_index]
                out_fname = "halos_mcut_{}_{}_{}_sample.txt".format(format_e(ml_cut), format_e(mu_cut), format_e(rand_sample))
                outfile_file.write(out_fname+'\n')
                print(out_fname)
        else:
                out_fname = "halos_mcut_{}_{}_.txt".format(format_e(ml_cut), format_e(mu_cut))
                outfile_file.write(out_fname+'\n')
                print(out_fname)

        outfile_file.close()

        data = data[(data['N'] > ml_cut) & (data['N'] <= mu_cut)]

        print(len(data))

        data.to_csv(output_path+out_fname, sep = ',', mode = 'a', index = False, header = False)
        return len(data)

#---------------------------------------------------------------------------

def load_files():
        # To do -> Create seperate mass list
        halo_files = os.listdir(file_path)
        for file in halo_files:
                if file[0:4] == 'halo':
                        halo_cat = readdata(file_path+file, ['N', 'SO_central_particle'])
                        halo_cat = halo_cat.sort_values('N')
                        halo_num = int(5e5) 
                        indeces = [i * halo_num for i in range(int(len(halo_cat)//halo_num)+1)]+[-1]
                        mass_list = halo_cat['N'].iloc[indeces].tolist()
                        for i in range(len(mass_list)-1):
                                create_mass_slice(halo_cat, mass_list[i], mass_list[i+1])

halo_cat = readdata('/home/vipul/vipul/halo_clutering/AbacusSummit_huge_c000_ph201/halos/z3.000/halo_info/halo_info_000.asdf', ['N', 'SO_central_particle']) # N = particle number
#halo_cat = halo_cat.sort_values('N')
'''
halo_num = int(3e5) # number of halos in each mass bin
indeces = [i * halo_num for i in range(int(len(halo_cat)//halo_num)+1)]+[-1]
mass_list = halo_cat['N'].iloc[indeces]
print(mass_list) #data_summary(halo_cat)

#for i in range(len(mass_list)-1):
#       create_mass_slice(halo_cat, mass_list[i], mass_list[i+1])
        
#m_bin_cat = create_mass_slice(halo_cat, 1e11, 5e11)
'''
