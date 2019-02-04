import numpy as np
import units_conversions as un
import parameters as par
import Cosmo_functions as cf
import BINGO_parameters as bpar
import Eisenstein_Hu_fits as EH
import HI_functions as hi
import scipy.integrate as integrate
from math import pi
from scipy.special import spherical_jn
from scipy.integrate import quad
from time import time

#==================================================================================
#==================================================================================
#= 3D Power Spectrum ==============================================================
#==================================================================================
#==================================================================================


def Pk_HI(k, z,params,info_spec):
	Dz     = cf.Dz_LCDM(z,params)
	bHI    = hi.bias(0,info_spec['bias'])
	Pk     = EH.Pk(k,params,info_spec['spectrum'])
	A      = EH.Amp_Pk(params,info_spec)
	z_type = un.verification_dtype_list(z)
	#print "Amplitude" + str(A/1.e6)

	if len(z_type)==1:return A*((bHI*Dz)**2)*Pk
	elif len(z_type)>1:
		Pkz = []
		for i in range(len(z_type)): Pkz.append(A*((bHI*Dz[i])**2)*Pk)
		return Pkz
		


#==================================================================================
#==================================================================================
#= Angular Power Spectrum =========================================================
#==================================================================================
#==================================================================================

#==================================================================================
#= Limber Approximation ===========================================================
#==================================================================================

def Cl_l_HI_Limber(l,z_bin, params, info_spec):
	
	k_l    = lambda x: np.sqrt(l*(l+1))/cf.comoving_distance(x,params)
	dcl_dz = lambda x: cf.window_function(x,z_bin,params,info_spec['window'])*hi.mean_brightness_temperature(x, params,info_spec['OmegaHI_model'])/cf.comoving_distance(x,params)
	dcl_dz2 = lambda x: dcl_dz(x)**2
	dcl_dz3 = lambda x: dcl_dz2(x)*cf.Hubble_z(x,params)*Pk_HI(k_l(x),x,params,info_spec)
	cl_l = integrate.quad(dcl_dz3,z_bin[0],z_bin[1])[0]/un.c_light
	return cl_l

def Cl_HI_Limber(l,z_bin, params, info_spec):
	cl = np.empty(len(l))
	for i,l_i in enumerate(l): 
		cl[i] = Cl_l_HI_Limber(l_i,z_bin, params, info_spec)
	return cl

def Dl_HI_Limber(l,z_bin, params, info_spec):
	return l*(l+1.)*Cl_HI_Limber(l,z_bin, params, info_spec)/2.*pi
