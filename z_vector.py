import numpy as np
import parameters as par
import BINGO_parameters as bpar





z_vector    = np.zeros(int(bpar.n_bins+1))
freq        = np.zeros(int(bpar.n_bins+1))
l_vector    = np.arange(2,400,1)
len_l       = len(l_vector)
DL          = np.zeros(len_l)
DLs         = np.zeros(len_l)
len_z       = len(z_vector)
z_vector[0] = bpar.z_min_BINGO


for n in range(len_z): z_vector[n] = (bpar.freq_21cm/(bpar.freq_min_BINGO + n*bpar.del_freq))-1.
z_vector = np.flip(z_vector,0)





def vec_names():
	for i in range(bpar.n_bins):
		if    i==0:names = ['bin1.dat']
		else: names.append('bin'+str(i+1)+'.dat')
	return names



