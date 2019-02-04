#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys,emcee
import numpy as np
import config_plots_cosmow as cpc
import MCMC_parameters as mpar
from matplotlib import pyplot as plt
from matplotlib import pylab, mlab, rc
from time import time
rc('text', usetex=True)
plt.rcParams['text.latex.unicode'] = True


########################################################################
#
# GetDist plots results
#
########################################################################

#print os.chdir(os.getcwd())
#print os.environ['HOME']
print "\n\n"
print "================================================================================================"
print "================================================================================================"
print "================================================================================================\n"
print "PROGRAMA: MCMC PARA PLOTAR OS RESAULTADOS DAS CADEIAS PARA OS PARAMETROS DO FIT\n"
print "================================================================================================"
print "================================================================================================"
print "================================================================================================"
print "\n\n"

########################################################################
#
# Type of plots
#
########################################################################

corner_plot  = False
getdist_plot = True



#========================================================================================================================
#=== Verificacao de pastas e cadeias ====================================================================================
#========================================================================================================================

#sys.exit(0)
if  len(os.listdir(os.path.join(os.getcwd(),'chains')))!=0:
	print "Existem cadeias para a(s) seguinte(s) divisao(oes) de bin(s):    " + str(os.listdir(os.path.join(os.getcwd(),'chains')))
	num_bins = int(raw_input("Qual divisao usar?  "))	
	print "\n"
	
else: 
	print"Nao ha cadeias. Primeiro, gere as cadeias"
	sys.exit(0)
	
print "Existem o(s) seguinte(s) bin(s):   " + str(len(os.listdir(os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins'))))

if  len(os.listdir(os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins')))!=0:
	bin_analysis = int(raw_input("Qual bin usar?  "))
	print "\n"
else: 
	print"Nao ha bins. Primeiro, gere o bin com as cadeias."
	sys.exit(0)	

num_chains = len(os.listdir(os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins','bin_'+str(bin_analysis))))-2
path_chains = os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins','bin_'+str(bin_analysis))
names_chains =  os.listdir(os.path.join(path_chains))
path_zeff    =  os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins','bin_'+str(bin_analysis),'zeff.dat')
path_chain0  =  os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins','bin_'+str(bin_analysis),'chain0.dat')


print "Em tal bin existem %i cadeias." %(num_chains)


#========================================================================================================================
#====  ==================================================================================================================
#========================================================================================================================
path_outputs = os.path.join(os.getcwd(),'outputs',str(num_bins)+'_bins')
if  os.path.exists(path_outputs)==False:
	os.mkdir(path_outputs)
else:
	if len(os.listdir(path_outputs))!=0:
		dirr = os.listdir(path_outputs)
		for filee in dirr:
			os.remove(path_chains_bin + os.sep + filee)
	else: pass
	


#========================================================================================================================
#========================================================================================================================
#========================================================================================================================
#for i in range(num_chains): print os.listdir(os.path.join(os.getcwd(),'chains',str(num_bins)+'_bins','bin_'+str(bin_analysis)))[i]

j=True
timei = time()
print "\n\n"
print "================================================================="
print " Carregando dados..."

for i in range(num_chains): 

	#if i>num_chains/3.:
		path_chain = os.path.join(path_chains,names_chains[i])
		if (path_chain!=path_chain0)*(path_chain!=path_zeff)==True:
			try:
				chain = np.loadtxt(path_chain)
				if j==True: 
					flat_sample = chain
					j=False
				else: flat_sample = np.vstack((flat_sample,chain))
	#			print path_chain, flat_sample.shape, chain.shape
					
			except: pass
			else: pass
print "\n"
print " Dados carregados."
print "================================================================="
print "\n\n"
	
del j,path_chain0,names_chains, num_chains

mylegends =  cpc.mylegends(np.loadtxt(path_zeff))

del path_zeff

print time()-timei
#sys.exit(0)
#========================================================================================================================
#========================================================================================================================
#========================================================================================================================


def print_datas_results(samples_f=None, chain_number_f = 1, getdist_sigma_f = cpc.getdist_sigma_error_bf):
	print "\n\n\n"
	print "================================================================="
	print "Resultados GetDist para a " + str(chain_number_f) + " cadeia :  "

	for j in range(len(cpc.mylabels)):
		print samples_f.getInlineLatex(cpc.mylabels[j],limit = int(getdist_sigma_f)) 
		#limit=1 (default) para regiao de 68%, limit=2 para a segunda regiao de contorno (especificada pelo  self.contours). VER: http://getdist.readthedocs.io/en/latest/mcsamples.html, GetInlineLatex().
	print "================================================================="
	print "\n\n"
'''
def save_datas(samples_f=None, chain_number_f = 1, getdist_sigma_f = cpc.getdist_sigma_error_bf,g_f = g, path_outputs_f = path_outputs):
	path_file_f = os.path.join(path_outputs_f,'mcmc_results.dat')
	
	print "\n\n\n"
	print "================================================================="
	print "Salvando os resultados e o plot em:  " + str(path_file_f)
	file_f = open(path_file_f,'a')
	for j in range(len(cpc.mylabels)):
		file_f.write(samples_f.getInlineLatex(cpc.mylabels[j]+ " " + limit = int(getdist_sigma_f)) + "\n")
	file_f.close()
	g_f.export(os.path.join(path_outputs,'triangle_plot_'+str(num_bins)+'.pdf'))
	print "Salvo"
	print "================================================================="
	print "\n\n"	
'''
#========================================================================================================================
#========================================================================================================================
#========================================================================================================================



if getdist_plot==True:	
	import getdist
	from getdist import plots, MCSamples
	
	g = plots.getSubplotPlotter()
	samples_list = [MCSamples(samples = flat_sample  ,names = cpc.mylabels , label = mylegends[0])]
	g.triangle_plot(samples_list, filled_compare = cpc.contour_filled, legend_labels = mylegends)
	g.settings.solid_contour_palefactor = 0.8
	g.settings.legend_fontsize = 10
	g.settings.alpha_filled_add= 0.4
	g.settings.figure_legend_frame = False
	
	if cpc.print_terminal_datas==True: print_datas_results(samples_list[0])
	#if cpc.save_datas==True:           save_datas(samples_f = samples_list[0], g_f = g, path_outputs_f = path_outputs)
	
	plt.show()
	

elif getdist_plot != False:
	print "\n\n"
	print "================================================================="
	print "Cuidado com as condiçôes de plotagem."
	print "================================================================="
	print "\n\n"
else: pass





########################################################################
#
# Corner plots results
#
########################################################################


if corner_plot==True:
	import corner
	print flat_sample.shape
	fig=corner.corner(flat_sample, labels = cpc.mylabels,truths=mpar.a_fid, show_titles=True, title_fmt=".3f")
	plt.show()
elif corner_plot != False:
	print "Cuidado com as condiçoes de plotagem."
else: pass
print time()-timei
