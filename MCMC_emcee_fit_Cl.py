#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os,sys,emcee
import numpy as np

import z_vector as z
import parameters as par
import BINGO_parameters as bpar
import MCMC_parameters as mpar
import Eisenstein_Hu_fits as EH
import Cosmo_functions as cf
from time import time

print "\n\n"
print "================================================================================================"
print "================================================================================================"
print "================================================================================================\n"
print "PROGRAMA: MCMC PARA CALCULAR OS PARAMETROS DO FIT FENOMENOLÓGICO DO Cl PARA OS BAOs\n"
print "================================================================================================"
print "================================================================================================"
print "================================================================================================"
print "\n\n"



########################################################################
#
# Way of inputs and outputs
#
########################################################################

nbins      = mpar.n_bins

bin_i      = mpar.bin_analysis
path       = os.path.join(os.getcwd(),'cl_bao',str(nbins) + '_bins')
path_bin_i = os.path.join(path,'bin'+str(int(bin_i))+'.dat')



########################################################################
#
# DATAS/INPUT
#
########################################################################


l,dl,dls = np.loadtxt(path_bin_i, unpack = True)
bao      = dl/dls
err      = np.sqrt(2./(2.*l+1.))*bao



print "================================================================="
print "Maior erro:   " + str(int(1e4*np.amax(err))/1e4)
print "Menor erro:   " + str(int(1e4*np.amin(err))/1e4)
print "================================================================="
print "\n\n"





########################################################################
#
# Define EMCEE parameters
#
########################################################################




print "================================================================="
print "Num. loops:   " + str(mpar.nsteps)
print "Num. passos:  " + str(mpar.steps)
print "Num. walkers: " + str(mpar.nwalkers)
print "Num. burning: " + str(mpar.burn_steps)
print "Num. pontos:  " + str(mpar.nwalkers*mpar.nsteps*mpar.steps)
print "================================================================="
print "\n\n"






########################################################################
#
# Vector of binned redshits values and sound horizon rs
#
########################################################################

z_vector  = np.loadtxt(os.path.join(path,'z.dat'))
zeff      = z_vector[bin_i]
rs        = EH.sound_horizon_drag(mpar.params)
#mylegends = cpc.mylegends(zeff)


print "================================================================="
print "Faixa de redshift do survey: (" + str(int(np.amin(z_vector)*1000.)/1000.) + "," + str(int(np.amax(z_vector)*1000.)/1000.) + ")"
print "Num. de bins:                 " + str(nbins)
print "Bin analisado:                " + str(bin_i)
print "Redshift do bin analisado:    " + str(int(zeff*1000.)/1000.)
print "================================================================="
print "\n\n"






########################################################################
#
# Define functions: Likelihood, Prior, data theory
#
########################################################################

def theory(a,l):

	A         = a[0]
	alpha     = a[1]
	SigmaSilk = a[2] 
	
	#kr   = l +0.5 +(-(1./(8.*l)) + (1./(16.*(l**2))) )
	kr   = np.sqrt(l*(l+1))
	k    =  kr/cf.comoving_distance(zeff,mpar.params)
	f_l  =  1. + A*k*np.exp(-(k*SigmaSilk)**1.4)*np.sin(alpha*k*rs)
	return f_l

def lnprior(a_try):
	if ((prior_range[:,0]<a_try)*(prior_range[:,1]>a_try)).all()==False: 
		return -np.inf
	return 0.0

def lnprob(a_try,data_vec):
    lp = lnprior(a_try) 
    if not np.isfinite(lp): 
        return -np.inf      
    try_vec = theory(a_try,l)
    diff_vec = try_vec - data_vec
    chi2_try = np.sum(diff_vec**2/err**2)	
    return -chi2_try/2.0

########################################################################
#
# 
#
########################################################################

def file_chains():
	path_file_bin   = os.path.join(os.getcwd() , 'chains' , str(nbins) + '_bins')
	path_chains_bin = os.path.join(path_file_bin,'bin_' + str(bin_i))
	
	if os.path.isdir(path_file_bin):
		if os.path.isdir(path_chains_bin):
			dirr = os.listdir(path_chains_bin)
			for filee in dirr:
				os.remove(path_chains_bin + os.sep + filee)
		else:
			os.mkdir(path_chains_bin)
			
	else:
		os.mkdir(path_file_bin)
		os.mkdir(path_chains_bin)

	g = open(os.path.join(path_chains_bin,'zeff.dat'),'w')
	g.write(str(zeff))
	g.close()

	return path_chains_bin



########################################################################
#
# Convergence Diagnostic: Gelman-Rubin
#
########################################################################

def Gelman_Rubin(chains_f=None):
	#calculate mean to mth chains:     Sigma_hat_m
	#calculate mean between chains:    Sigma_hat
	#calculate variance to mth chains: Var_hat_m
	#Let N be elements chains number in a chains
	#Let M be chains number
	#B = frac(N)(M-1) Sum_{m=1}^{N} (Sigma_hat_m - Sigma_hat)**2 , between-chains variance
	#W = frac(1)(M)Sum_{m=1}^{M} (Var_hat_m)^{2} ,                 within-chains variance
	#V_hat = frac(N-1)(M)W + frac(M+1)(NM)B
	return None

########################################################################
#
# Initial chains and conditions
#
########################################################################

prior_range = np.array([mpar.a0_range,mpar.a1_range,mpar.a2_range])


a0_0  = np.random.uniform(mpar.a0_range[0],mpar.a0_range[1],mpar.nwalkers)
a1_0  = np.random.uniform(mpar.a1_range[0],mpar.a1_range[1],mpar.nwalkers)
a2_0  = np.random.uniform(mpar.a2_range[0],mpar.a2_range[1],mpar.nwalkers)
par_0 = np.array([a0_0,a1_0,a2_0]).T 


########################################################################
#
# Begin program
#
########################################################################

sampler = emcee.EnsembleSampler(mpar.nwalkers, mpar.ndim, lnprob, args=[bao])

timei=time()
pos, prob, state = sampler.run_mcmc(par_0, mpar.burn_steps)  
burnsamp = sampler.chain 
burn = range(mpar.burn_steps) 
timef=time()
delta_t=(timef-timei)/mpar.nwalkers/mpar.burn_steps 

print "================================================================="
print("Time for each MCMC point:", delta_t)


path_chains = file_chains()
f = open(os.path.join(path_chains,"chain0.dat"), "w")
f.close()

sampler.reset() 
posf, probf, statef = sampler.run_mcmc(pos, mpar.steps)
flat_sample = sampler.chain.reshape((-1, mpar.ndim))




f = open(os.path.join(path_chains,"chain0.dat"), "a")
for k in range(flat_sample.shape[0]):
    f.write(str(flat_sample[k])[2:-2] + "\n")
f.close()


print("Estimated duration of MCMC:", delta_t*mpar.nwalkers*mpar.nsteps*mpar.steps)
print "================================================================="
print "\n\n"

print "================================================================="
print ("Criando as cadeias...")

for nstore in range(mpar.nsteps):
    
    posf, probf, statef = sampler.run_mcmc(posf, mpar.steps)
    f = open(os.path.join(path_chains,"chain"+str(nstore+1)+ ".dat"), "a")
    this_sample = sampler.chain.reshape((-1, mpar.ndim))
    for k in range(this_sample.shape[0]):
        f.write(str(this_sample[k])[2:-2] + "\n")
    f.close()
#print("\n")
print("Cadeias criadas.")
print "================================================================="
