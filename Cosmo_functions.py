#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.integrate as integrate
import parameters as par
import units_conversions as un
from math import pow,sqrt,pi,sin,sinh
from scipy.special import spherical_jn


#=======================================================================
#=======================================================================
# Equation of states ===================================================
#=======================================================================
#=======================================================================

def w_function(z,params):
	if params['model_w']=='const':
		return params['w']
		
	elif params['model_w']=='cpl':
		w = params['wa']
		w = w*z
		w = w/(1. + z)
		return paras['w'] + w
	
	elif params['model_w']=='Interaction':
		print('No implemented')
		sys.exit(0)
		return None

	elif params['model_w']=='quintessence':
		print('No implemented')
		sys.exit(0)
		return None
				
	else:
		print('Erro')
		print('End program !')
		sys.exit(0)
		
#=======================================================================
#=======================================================================
# Hubble's functions H(z) ==============================================
#=======================================================================
#=======================================================================

def Hubble_a(a, params): 
	
	h2 = params['h']**2
	w  = w_function((-1 + 1./a),params)
	Ok = params['Ok_h2']/h2/pow(a,2)
	Om = params['Om_h2']/h2/pow(a,3)
	Or = params['Or_h2']/h2/pow(a,4)
	Od = params['Od_h2']/h2/pow(a,-3.*(1+w)) 
		
	return 100*params['h']*sqrt(Om + Or + Ok + Od)

def Hubble_z(z, params): 
	
	h2 = params['h']**2
	w  = w_function(z,params)
	Ok = params['Ok_h2']/h2*pow(1+z,2)
	Om = params['Om_h2']/h2*pow(1+z,3)
	Or = params['Or_h2']/h2*pow(1+z,4)
	Od = params['Od_h2']/h2*pow(1+z,3.*(1+w))
	return 100*params['h']*sqrt(Od + Ok + Om + Or)

#=======================================================================
#=======================================================================
# Distance functions  ==================================================
#=======================================================================
#=======================================================================

def luminosity_distance(z, params):
	
	H0       = 100*params['h']
	integral = lambda x: H0/Hubble_z(x,params)
	integral = integrate.quad(integral,0.,z)[0]
	
	if params['Ok_h2']==0:
		Dl = un.c_light*(1.+z)
		Dl = Dl*integral/H0
		return Dl
		
	elif params['Ok_h2']<0:
		h2 = params['h']**2
		Dl = sqrt(-params['Ok_h2']/h2)*integral
		Dl = sin(Dl)
		Dl = Dl/sqrt(-params['Ok_h2'])/100.
		Dl *=un.c_light*(1+z)
		return Dl		
	elif params['Ok_h2']>0:
		h2 = params['h']**2
		Dl = sqrt(params['Ok_h2']/h2)*integral
		Dl = sinh(Dl)
		Dl = Dl/sqrt(params['Ok_h2'])/100.
		Dl *=un.c_light*(1+z)
		return Dl
			
	else: 
		print('error')
		import sys
		sys.eixt(0)


def angular_distance(z, params):
	da = luminosity_distance(z,params)
	da = da/((1.+z)**2)
	return da


def comoving_distance(z, params):
	dc = lambda x: 1./Hubble_z(x,params)
	dc = integrate.quad(dc,0.,z)[0]
	dc = un.c_light*dc
	return dc 

def angle_averaged_distance(z, params):
	H   = Hubble_z(z, params)
	da  = angular_distance(z, params)
	da  = da*(1.+z)
	da  = da**2
	dv  = da*z
	dv  = un.c_light*dv/H
	dv  = pow(dv,1./3.)
	return dv

#=======================================================================
#=======================================================================
#= Growth Function =====================================================
#=======================================================================
#=======================================================================

def Dz_z_LCDM(z, params):#Dodelson, Modern Cosmology pg 207-211, eq. 7.77, for to no interation models
	
	if params['Ok_h2']==0 or params['Om_k']==0.: # Growth_function. This function is to flat universe.
		a_z      = 1./(1. + z)
		integral = integrate.quad(lambda x: 1./pow(x*Hubble_a(x,params)/(100*params['h']),3),0.,a_z)[0]
		Om       = params['Om_h2']/(params['h']**2)
		H_a      = Hubble_z(z,params)/(100*params['h'])
		
		D1       = 5.*Om*H_a/2.
		D1       = D1*integral
		
		norm     = integrate.quad(lambda x: 1./pow(x*Hubble_a(x,params)/(100*params['h']),3),0.,1.)[0]
		norm     = 5.*Om*norm/2.
		
		return D1/norm
		
	elif params['Ok_h2']>0: # Growth_function. This function is to open universe. k<0 -> Ok>0
		a_z = 1./(1. + z)
		Om  = params['Om_h2']/(params['h']**2)
		x   = 1. - Om
		x   = x*a_z/Om
		
		D1  = log(sqrt(1+x) - sqrt(x))
		D1  = 3.*D1*sqrt(1+x)
		D1  = D1/pow(x,3./2.)
		D1  += 3./x
		D1  += 1
		D1  *= (5.*a_z)/(2.*x)
		
		x    = (1. - Om)/Om
		norm = log(sqrt(1+x) - sqrt(x))
		norm = 3.*D1*sqrt(1+x)
		norm = D1/pow(x,3./2.)
		norm += 3./x
		norm += 1
		norm *= 5./(2.*x)
		
		return D1/norm
	else:
		print("Erro")
		print("Finish program !")
		sys.exit(0)
	return None

def Dz_LCDM(z,params):
	z_type = un.verification_dtype_list(z)
	if len(z_type)>1: 
		D = np.empty(len(z))
		for i,z_i in enumerate(z): D[i] = Dz_z_LCDM(z_i,params)
		return D
	elif len(z_type)==1: return Dz_z_LCDM(z,params)
	else:
		print("Erro")
		print("Finish program")
		sys.exit(0)
		

		
#=======================================================================
#=======================================================================
#= Window Function =====================================================
#=======================================================================
#=======================================================================

def window_function(z,z_bin, params, window): #window = "battye","tophat","gaussian","zero","one"
	z_min = z_bin[0]
	z_max = z_bin[1]

	if window == "battye":
		if z_min==z_max:
			return 0
		elif z>=z_min and z<=z_max:
			W = z_max - z_min
			return 1./W
		else: 
			return 0.
	
	elif window =="tophat": 
		if z_min==z_max:# exijo que exista intervalo, senão não faz sentido o Intesity Mapping(IM)
			return 0
		elif z>=z_min and z<=z_max:
			return 1.
		else:
			return 0.
	
	elif window == "gaussian":
		kR  = k*params['R']
		kR2 = kR**2
		WR  = exp(-kR2/2.)
		WR  = WR/((sqrt(2.*pi)*params['R'])**3)		
		return WR
			
	elif window == "zero":
		return 0.
	
	elif window == "one":
		return 1.
	
	else:
		print("Erro")
		print("Finish Program")
		import sys
		sys.exit(0)
		
	return None

def fourier_window_function(k, params, symmetry, window):
	if symmetry == "spherical" and window == "tophat":
		kR  = k*params['R']
		kR3 = kR**3
		WR  = (kR)*cos(kR)
		WR  = sin(kR) - WR
		WR  = 3.*WR
		return WR/(kR3)
		
	elif symmetry == "spherical" and window == "gaussian":
		kR  = k*params['R']
		kR2 = kR**2
		WR  = exp(-kR2/2.)
		return WR
	
	else:
		print("Erro")
		print("Finish Program")
		sys.exit(0)
	
	return None

#=======================================================================
#=======================================================================
#= Parameters ==========================================================
#=======================================================================
#=======================================================================

def AP_parameter(z,params): #Alcock-Paczynski parameter
	F_ap = (1+z)/un.c_light
	F_ap *= angular_distance(z,params)*Hubble_z(z, params)
	return F_ap
	
def acoustic_parameter_A(z, params): #Eisenstein et al. 2005
	A = sqrt(params['Om_h2']*(1.e4))/un.c_light
	A *= angle_averaged_distance(z, params)/z
	return A
	
