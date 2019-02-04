import BINGO_parameters as bpar
import numpy as np
c_light_m_s        = 2.99792458e8          #m/s
c_light            = c_light_m_s/1e3       #km/s
G_const_Mpc_Msun_s = 4.51737014558*1.e-48
G_const            = 6.67408*1.e-11        #(m3)/(kg*s2)
G_const_km         = 6.67408*1.e-2         #(km3)/(kg*s2)
Mpc_cm             = 3.08568025*1.e+24     #Mpc in cm
Mpc_km             = 3.08568025*1.e+19     #Mpc in km
K4_gcm3            = 1.279*10**(-35)       # 1 k**4 = 1.279.1e-35 g cm-3

def convers_z_nu(z=None, nu_emit = bpar.freq_21cm):
	return nu_emit/(1.+z)
	
def convers_nu_z(nu=None, nu_emit = bpar.freq_21cm):
	return (nu_emit/nu)-1.

def convers_lamb_nu(lamb = bpar.Lambda_21cm):
	return par.c_light*1.e3/lamb

def convers_nu_lamb(nu = bpar.freq_21cm):
	return par.c_light*1.e3/nu

def verification_dtype_list(x):
	
	if (type(x)==int) or (type(x)==float):
		return[x]
	elif type(x)==list:
		return x
	elif (x.dtype == np.dtype(float))or(x.dtype == np.dtype(int)):
		x = x.tolist()
		if (type(x)==float)or(type(x)==int):
			return [x]
		else:
			return x
	else: 
		return x	


