

import numpy as np
import parameters as par
import scipy.integrate as integrate
import units_conversions as un
from math import *#log,e,sqrt,pow
from scipy.special import spherical_jn


#==================================================================================
#==================================================================================
#= Cosmological variables =========================================================
#==================================================================================
#==================================================================================

def z_eq(params):
	zeq = params['Om_h2']/params['Or_h2']
	zeq = zeq-1.
	return zeq

def k_eq(params): 
	keq = 2.*params['Om_h2']*(100**2)
	keq = keq*z_eq(params)
	keq = sqrt(keq)
	return keq/un.c_light

def k_silk(params):
	ksilk   = pow(10.4*params['Om_h2'],-0.95)
	ksilk   = 1.+ ksilk
	ksilk   = pow(params['Om_h2'],0.73)*ksilk
	ksilk   = pow(params['Ob_h2'],0.52)*ksilk
	return 1.6*ksilk #Mpc-1

#==================================================================================
#==================================================================================
# Fits variables for build of the functions transfer ==============================
#==================================================================================
#==================================================================================

def b1(params):
	b_1 = 0.607*pow(params['Om_h2'],0.674)
	b_1 = 1. + b_1
	b_1 = b_1*0.313/pow(params['Om_h2'],0.419)
	return b_1

def b2(params):
	return 0.238*pow(params['Om_h2'],0.223)

def z_drag(params):
	b_1 = b1(params)
	b_2 = b2(params)
	zd = 1291*pow(params['Om_h2'],0.251)
	zd = zd*(1. + b_1*pow(params['Ob_h2'],b_2))
	zd = zd/(1. + 0.659*pow(params['Om_h2'],0.828))
	return zd

def R_parameter(z,params):
	R = (3.*params['Ob_h2'])/(4.*params['Og_h2'])
	R = R/(1.+z)
	return R

def sound_speed(z,params): 
	R = R_parameter(z,params)
	return un.c_light/sqrt(3.*(1.+R)) #km.s-1

def sound_horizon_drag(params): 
	
	zeq     = z_eq(params)
	zdrag   = z_drag(params)
	keq     = k_eq(params)
	R_drag  = R_parameter(zdrag,params)
	R_eq    = R_parameter(zeq,params)
		
	s       = log((sqrt(1.+R_drag) + sqrt(R_drag + R_eq))/(1. +sqrt(R_eq)))
	s       = sqrt(6./R_eq)*s
	s       = (2./(3.*keq))*s
	return s 


#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================

def s_parameter	(params):
	s = 44.5*log(9.83/params['Om_h2'])
	s = s/sqrt(1. + 10.*pow(params['Ob_h2'],0.75))
	return s

def k_peak(params):
	s = s_parameter(params)
	k = 5*pi*(1. + 0.217*params['Om_h2'])
	return k/(2.*s)

def alpha_gamma(params):
	f     = params['Ob_h2']/params['Om_h2']
	f2    = f**2
	alpha = 0.38*log(22.3*params['Om_h2'])*f2
	alpha = alpha - 0.328*log(431.*params['Om_h2'])*f 
	return 1. + alpha

def gamma(params):
	return params['Om_h2']/params['h']
	
def gamma_eff(ks, params):
	alpha  = alpha_gamma(params)
	return gamma(params)*(alpha + ((1.- alpha)/(1. + (0.43*ks)**4)))

def q_parameter(k, gamma, params):
	theta  = params['TCMB']/2.7
	theta2 = theta**2
	return (k/params['h'])*(theta2/gamma) #k should are in h Mpc-1.

def C0(q):
	c = 731./(1.+62.5*q)
	return c+14.2

def L0(q): 
	return log(2.*e +1.8*q)

def T0(q):
	t=L0(q)
	t=t/(t+(C0(q)*(q**2)))
	return t

#==================================================================================
#==================================================================================
# Fit Eisenstein&Hu ===============================================================
#==================================================================================
#==================================================================================

def a_1(params):
	a       = 32.1*params['Om_h2']
	a       = 1.+ pow(a,-0.532)
	a       = a*pow(46.9*params['Om_h2'],0.670)
	return a

def a_2(params):
	a       = 45.0*params['Om_h2']
	a       = 1.+ pow(a,-0.582)
	a       = a*pow(12.0*params['Om_h2'],0.424)
	return a
	
def b_1(params):
	b       = 458.*params['Om_h2']
	b       = 1.+pow(b,-0.708)
	b       = 0.944/b
	return b

def b_2(params):
	b       = 0.395*params['Om_h2']
	b       = pow(b,-0.0266)
	return b

def alpha_b(params):
		
	zeq     = z_eq(params)
	zd      = z_drag(params)
	Rd      = R_parameter(zd,params)
	keq     = k_eq(params)
	s       = s_parameter(params)
	ks      = keq*s
	G       = G_function(zeq,zd)
	
	alpha = 2.07*ks*G
	alpha = alpha/pow(1.+Rd,0.75) 
	return alpha

def alpha_c(params):
	f      = params['Ob_h2']/params['Om_h2']
	f3     = f**3
	a1      = pow(a_1(params),-f)
	a2      = pow(a_2(params),-f3)
	return a1*a2

def beta_b(params):
	f       = params['Ob_h2']/params['Om_h2']
	beta    = 1.+ (17.2*params['Om_h2'])**2
	beta    = sqrt(beta)
	beta    = (3. - 2.*f)*beta
	beta    = 0.5 + f + beta
	return beta

def beta_c(params):
	fc      = params['Oc_h2']/params['Om_h2']
	b1      = b_1(params)
	b2      = b_2(params)
	beta    = (pow(fc,b2)-1.)
	beta    = b1*beta
	beta    = beta +1.
	return 1./beta

def beta_node(params):
	return 8.41*pow(params['Om_h2'],0.435)

def G_function(zeq, zd):
	x = (1.+zeq)/(1.+zd)
	G = (sqrt(1.+x)+1.)/(sqrt(1.+x)-1.)
	G = (2.+3.*x)*log(G)
	G = G - 6.*sqrt(1.+x)
	G = x*G
	return G

def T0_til(q, k, alphac, betac):
	q2 = q**2
	C  = (1.+69.9*pow(q,1.08))
	C  = 386./C
	C  = C+(14.2/alphac)
	
	T  = log(e + 1.8*betac*q) 
	T  = T/(T+C*q2)
	return T

def s_shifting(k,params):
	if k==0 or k==0.:return 0.
	betan      = beta_node(params)
	s          = s_parameter(params)
	ks         = k*s
	betan      = (betan/ks)**3
	s_til      = pow(1.+ betan,1./3.)
	return s/s_til

#==================================================================================
#==================================================================================
#= Transfer Functions : Baryons and CDM ===========================================
#==================================================================================
#==================================================================================

def Tb_eisenstein(k,params):
	if k==0 or k==0.: return 1.
	
	alphab   = alpha_b(params)
	betab    = beta_b(params)
	
	q        = q_parameter(k, gamma(params), params)
	ksilk    = k_silk(params)
		
	s        = s_parameter(params)
	s_til    = s_shifting(k, params)
	ks       = k*s
		
	T0       = T0_til(q,k,1.,1.)
	Tb       = T0/(1.+(ks/5.2)**2)
	Tb       = Tb + (alphab/(1.+pow(betab/ks,3)))*exp(-pow(k/ksilk,1.4))
	
	j0       = spherical_jn(0,k*s_til,0)
	Tb       = Tb*j0
	return Tb

def Tc_eisenstein(k, params):
	ks      = k*s_parameter(params)
		
	alphac  = alpha_c(params)
	betac   = beta_c(params)
	q       = q_parameter(k, gamma(params),params)
	
	T0      = T0_til(q,k,alphac,betac)
	T0_1    = T0_til(q,k,1.,betac)
	f       = 1. + pow(ks/5.4,4)
	f       = 1./f
	
	Tcdm    = f*T0_1 + (1.-f)*T0
	return Tcdm

#==================================================================================
#==================================================================================
#= Transfer Functions : NoBAO, Nobaryons, BBKS and Eisenstein&HU ==============
#==================================================================================
#==================================================================================

def Tm_nobao(k, params):
	ks    = k*s_parameter(params)
	Gamma = gamma_eff(ks, params)
	q_eff = q_parameter(k, Gamma, params)
	return T0(q_eff)

def Tm_nobaryon(k, params):
	q = q_parameter(k, gamma(params), params)
	return T0(q)

def Tm_BBKS(k, params):# fit Bardeen,Bond,Kaiser & Szalay(BBKS)
	keq = k_eq(params)
	x   = k/keq
	TF  = 1.+(0.284*x)+pow(1.18*x,2)+pow(0.399*x,3)+pow(0.490*x,4)
	TF  = pow(TF,-0.25)
	TF  = TF*log(1.+0.171*x)
	if k!=0:TF = TF/(0.171*x)
	else: TF=1.
	return TF

def Tm_eisenstein(k, params):
	                                                                                               
	fb = params["Ob_h2"]/params["Om_h2"]
	fc = params["Oc_h2"]/params["Om_h2"]
	Tc = Tc_eisenstein(k, params)
	Tb = Tb_eisenstein(k, params)
	return fb*Tb + fc*Tc 

#==================================================================================
#==================================================================================
#= Power Spectrum : NoBAO, Nobaryons, BBKS and Eisenstein&HU ==============
#==================================================================================
#==================================================================================

def Pk(k,params,spectrum): # without ampl #spectrum = "eisenstein","nobao","nobaryon","bbks"

	k_type = un.verification_dtype_list(k)
	if len(k_type)>1:
		TF = np.empty(len(k))
		
		if spectrum == "eisenstein":
			for i,k_i in enumerate(k):
				TF[i] = Tm_eisenstein(k_i,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
			
		elif spectrum == "nobao":
			for i,k_i in enumerate(k):
				TF[i] = Tm_nobao(k_i,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
			
		elif spectrum == "nobaryon":
			for i,k_i in enumerate(k):
				TF[i] = Tm_nobaryon(k_i,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
			
		elif spectrum == "bbks":
			for i,k_i in enumerate(k):
				TF[i] = Tm_BBKS(k_i,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
		else: 
			print("Error")
			print("Finish program !!!")
			sys.exit(0)
	
	
	elif len(k_type)==1:
		
		if spectrum == "eisenstein":
			TF = Tm_eisenstein(k,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
			
		elif spectrum == "nobao":
			TF = Tm_nobao(k,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
						
		elif spectrum == "nobaryon":
			TF = Tm_nobaryon(k,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
		
			
		elif spectrum == "bbks":
			TF = Tm_BBKS(k,params)
			p_ini = k**params['ns']
			return p_ini*(TF**2)
		
		else: 
			print("Error")
			print("Finish program !!!")
			sys.exit(0)
		
def dimensionless_Pk(k,params,spectrum): #without ampl.
	TF = np.empty(len(k))
	if spectrum == "eisenstein":
		for i,k_i in enumerate(k):
			TF[i] = Tm_eisenstein(k_i,params)
		p_ini = k**(params['ns']+3)
		return p_ini*(TF**2)/(2.*pi**2)
		
	elif spectrum == "nobao":
		for i,k_i in enumerate(k):
			TF[i] = Tm_nobao(k_i,params)
		p_ini = k**(params['ns']+3)
		return p_ini*(TF**2)/(2.*pi**2)
		
	elif spectrum == "nobaryon":
		for i,k_i in enumerate(k):
			TF[i] = Tm_nobaryon(k_i,params)
		p_ini = k**(params['ns']+3)
		return p_ini*(TF**2)/(2.*pi**2)
		
	elif spectrum == "bbks":
		for i,k_i in enumerate(k):
			TF[i] = Tm_BBKS(k_i,params)
		p_ini = k**(params['ns']+3)
		return p_ini*(TF**2)/(2.*pi**2)
		
	else: 
		print("Error")
		print("Finish program !!!")
		sys.exit(0)

#=======================================================================
#=======================================================================
#= Amplitude with sigmaR ===============================================
#=======================================================================
#=======================================================================
def delta_H(params):
	n_til = params['ns'] - 1.
	h     = params['h']
	Om    = params['Om_h2']/h**2
	if params['Ok_h2']==0 or params['Ok_h2']==0.:
		delH = np.exp(-(n_til + 0.14*n_til**2))
		delH *= Om**(-(0.35 + 0.19*log(Om)+0.017*n_til))
		delH *= 1.95e-5
		return delH
	elif params['Ok_h2']>0:
		delH = np.exp(-(0.95*n_til + 0.169*n_til**2))
		delH *= Om**(-(0.785 + 0.05*log(Om)))
		delH *= 1.94e-5
		return delH
	else: 
		print("Error")
		print("Finish program !!!")
		sys.exit(0)


def Amp_Pk(params, info_spec):
	
	if info_spec['amplitude']=='sigma8':
		dPkj1_dk = lambda x: Tm_eisenstein(x,params)*(spherical_jn(1,x*params['R'],0)**2)
		A = integrate.quad(dPkj1_dk,0,np.inf)[0]
		A *= (9./2./(pi*params["R"])**2)
		return params["sigma_R"]**2/A
	
	elif info_spec['amplitude']=='cobe':
		A = 2.*(pi*delta_H(params))**2
		A = A/(100*params['h']/un.c_light)**(3 + params['ns'])
		return A
		
	else: 
		print("Error")
		print("Finish program !!!")
		sys.exit(0)
