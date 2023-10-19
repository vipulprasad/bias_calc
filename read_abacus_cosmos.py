"""
Here halo mass is provided as number of particles in a halo,
position is given as -550 to +550 Mpc

 ParticleMassMsun =     5.411408436271152e+10
 ParticleMassHMsun =     3.731168281372334e+10
 
 BoxSizeMpc = 1595.358030249148442
 BoxSizeHMpc = 1100.000000000000000"""

import numpy as np
import random
import pandas as pd

#halo mass is provided as number of particles
particle_mass_real = 5.776643562630781e+10 #Msun
particle_mass =     3.731168281372334e+10 #hMsun
boxsize = 1100.000000000000000 #hMPc
boxsize_real = 1635.444543562294939 #MPc


#struct form of the stored data
halo_dt_FoF = lambda: np.dtype([("id",np.int64),("subsamp_start",np.uint64),("subsamp_len",np.uint32),
               ("N",np.uint32),("subhalo_N",np.uint32,4),("pos",np.float32,3),("vel",np.float32,3),("sigma_v",np.float32,3),
               ("r25",np.float32),("r50",np.float32),("r75",np.float32),("r90",np.float32),
               ("vcirc_max",np.float32),("rvcirc_max",np.float32),
               ("subhalo_pos",np.float32,3),("subhalo_vel",np.float32,3),("subhalo_sigma_v",np.float32,3),
               ("subhalo_r25",np.float32),("subhalo_r50",np.float32),("subhalo_r75",np.float32),("subhalo_r90",np.float32),
               ("subhalo_vcirc_max",np.float32),("subhalo_rvcirc_max",np.float32)], align=True)

def format_e(n):
    a = '%.2E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


def read_file(filename, rand = 1):
    # rand should be a value between 0 and 1, fraction for random sampling
	with open(filename,"rb") as fp:
		
		#reading to structs
		num_groups = np.fromfile(fp,dtype=np.uint64,count=1)
		n_largest_subhalos = np.fromfile(fp,dtype=np.uint64,count=1)
		assert (n_largest_subhalos,) == halo_dt_FoF()['subhalo_N'].shape, n_largest_subhalos
		data = np.fromfile(fp, dtype=halo_dt_FoF())
		assert len(data) == num_groups, (len(data), num_groups, halo_fn)
		
		#taking number of halo particles and position
		data = data[["N","pos"]]#N is the number of particles in a halos 
		
		if rand != 1:#ramodom sampling
		    random_data = []
		    N = len(data)
		    index = random.sample(range(N), int(rand*N)) # 0.1 gives 10% data
		    for i in index:
			    random_data.append(data[i])
		    data = random_data

		mass = []
		pos = []
		
		for i in range(len(data)):
			mass.append(data[i][0]*particle_mass)
			pos.append(data[i][1]+550) # +550 is to make the box coordinates between 0 and 1100Mpc
		
		del data
		return pd.DataFrame(np.column_stack((mass, pos)), columns = ['M', 'X', 'Y', 'Z'])	

def create_mass_slice(data, ml_cut, mu_cut):

        out_fname = "halos_mcut_{}_{}_.txt".format(format_e(ml_cut), format_e(mu_cut))
        print(out_fname)

        data = data[(data['M'] > ml_cut) & (data['M'] <= mu_cut)]

        print(len(data))

        data.to_csv(output_path+out_fname, sep = ',', mode = 'a', index = False, header = False)
        print(len(data))

input_path = '/home/vipul/vipul/halo_clutering/bias_calc/abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z0.300/'
output_path = '/home/vipul/vipul/halo_clutering/bias_calc/abacus_cosmos/AbacusCosmos_1100box_planck_00-0_FoF_halos/z0.300/'
mass_bins = [1e12, 1.5e12, 2e12, 3e13, 2e15]

for i in range(4):
	file_name = input_path+'halos_{}'.format(i)
	halo_cat = read_file(file_name)
	create_mass_slice(halo_cat, 1e12, 1.5e12)
	create_mass_slice(halo_cat, 2e12, 3e13)
	create_mass_slice(halo_cat, 3e13, 2e15)

