import units_conversions as un
import parameters as par
import Cosmo_functions as cf
import BINGO_parameters as bpar
from math import pi


def f_sky_BINGO(V_bingo, unit): #V_bingo=bpar.Omega_Survey, unit = "deg"
	
	if unit == "deg":V_sky = 4.*pi*((180./pi)**2) #f_sky = 4pi**2 [rad2] = 4pi(180/pi)**2 [deg2]
	else: V_sky = 4.*(pi**2)
	
	return V_bingo/V_sky
	

def l_max(nu_obs): return 2.*pi*bpar.D_antenna_diam/(un.c_light*10**3/nu_obs) #see: A.Hall et al, (2013) e R. Battye (2012)
                                                                              #see too: M.ZALDARRIAGA, S.FURLANETTO & L.HERNQUIST (2004)- arxiv: astro-ph/0311514v2

def Omega_pix(unit):
	if unit == "deg": return (bpar.Theta_FWHM/60.)**2 #deg2
	elif unit == "rad": return ((pi/180.)**2)*(bpar.Theta_FWHM/60.)**2 #rad2
	else: 
		print("Error")
		print("Finish program")
		import sys
		sys.exit(0)
		return None
	 
	


def t_pix(unit):     
	year_sec = 365.*24.*60.*60.
	tobs     = bpar.t_obs*year_sec #tempo em s.
	return bpar.n_f*bpar.tobs*Omega_pix(unit)/bpar.Omega_Survey	
