import sys,os
import numpy as np
#import units_conversions as un
import parameters as par
#import Cosmo_functions as cf
#import BINGO_parameters as bpar
#import Eisenstein_Hu_fits as EH
#import HI_functions as hi
import Limber_approx as cl
import matplotlib.pyplot as plt
import camb
from camb import model, initialpower
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['text.latex.unicode'] = True


params = {"R":par.R_sigma,"sigma_R":par.sigma8,"Om_h2":par.Om_h2,"Ob_h2":par.Ob_h2,"Oc_h2":par.Oc_h2,"Or_h2":par.Or_h2,"Og_h2":par.Og_h2, "Ok_h2":par.Ok_h2, "Od_h2":par.Od_h2,"h":par.h,"TCMB":par.TCMB,"ns":par.ns,'w':par.w,'model_w':par.model_w}
info_spec = {"spectrum": "eisenstein","amplitude": "cobe", "bias":"constant","OmegaHI_model":"constant", "info plot": "BAO","window":"battye"} #info_plot = 'BAO'#'single',smooth'
k  = np.logspace(-4,np.log10(3), num=200)
z=0
Pk =  cl.Pk_HI(k*params['h'],z,params,info_spec)

pars = camb.CAMBparams()

pars.set_cosmology(H0 = 100*params['h'], ombh2 = params['Ob_h2'], omch2 = params['Oc_h2'], omk = params['Ok_h2']/params['h']**2)
pars.InitPower.set_params(ns=params['ns'],r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0)
pars.set_dark_energy()
pars.set_matter_power(redshifts=[0], kmax=3.0)
results = camb.get_results(pars)
#data = camb.get_transfer_functions(pars)
#transfer = data.get_cmb_transfer_data()

kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
plt.loglog(kh,pk[0,:],color='blue', label='camb')
A = pk[0][0]/Pk[0]
plt.loglog(k,A*Pk, color = "red",label = "Eisenstein\&Hu")
plt.xlim(kh[0],np.amax(kh))
plt.xlabel(r'$k\,\left(Mpc^{-1}h^{-1}\right)$')
plt.ylabel(r"$P\left(k\right)\,\left(Mpc/h\right)^3$")
plt.legend()
plt.show()




