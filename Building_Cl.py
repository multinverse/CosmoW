import os, sys
import os.path
import commands
import numpy as np
import z_vector as z
import parameters as par
import BINGO_parameters as bpar
import Limber_approx as li 
from time import time



print "\n\n"
print "================================================================================================"
print "================================================================================================"
print "================================================================================================\n"
print "PROGRAMA: GERAR ESPECTROS DE Cl e Cls DE 21cm E SALVA-LOS EM ARQUIVOS EXTERNOS\n"
print "================================================================================================"
print "================================================================================================"
print "================================================================================================"
print "\n\n"



params = {"R":par.R_sigma,"sigma_R":par.sigma8,"Om_h2":par.Om_h2,"Ob_h2":par.Ob_h2,"Oc_h2":par.Oc_h2,"Or_h2":par.Or_h2,
          "Og_h2":par.Og_h2, "Ok_h2":par.Ok_h2, "Od_h2":par.Od_h2,"h":par.h,"TCMB":par.TCMB,"ns":par.ns,'w':par.w,'model_w':par.model_w}
info_spec = {"spectrum": "eisenstein","amplitude": "cobe", "bias":"constant","OmegaHI_model":"constant","window":"battye"}


path = os.getcwd()
path_clbao = os.path.join(path , 'cl_bao')
path_output = os.path.join(path , 'output')
path_chains = os.path.join(path , 'chains')
path_bins  = os.path.join(path , 'cl_bao' , str(bpar.n_bins) + '_bins')


zz     = z.z_vector
names  = z.vec_names()
l      = np.arange(par.lmin,par.lmax+1,1)
l_save = l.tolist()


if os.path.isdir(path_clbao):
	pass
else:
	os.mkdir(path_clbao)

if os.path.isdir(path_output):
	pass
else:
	os.mkdir(path_output)
	
if os.path.isdir(path_chains):
	pass
else:
	os.mkdir(path_chains)
#=======================================================================
#=======================================================================
# verificar se a pasta n_bins, existe, e se ha arquivos nela
#=======================================================================
#=======================================================================

if os.path.isdir(path_bins):
	dirr = os.listdir(path_bins)
	for filee in dirr:
		os.remove(path_bins + os.sep + filee) 	
else:
	os.mkdir(path_bins)
	
#=======================================================================
#=======================================================================
#=======================================================================
#=======================================================================

print "Dados serao armazenados em:  " + str(path_bins)
print '\n\n'
for i in range(bpar.n_bins):
	t = time()
	info_spec['spectrum']='eisenstein'
	dl = li.Dl_HI_Limber(l,[zz[i],zz[i+1]],params, info_spec)
	info_spec['spectrum']='nobao'
	dls = li.Dl_HI_Limber(l,[zz[i],zz[i+1]],params, info_spec)
	
	name_file = path_bins + os.sep + names[i]
	
	#dl_save  = dl.tolist()
	#dls_save = dls.tolist()	
	#np.savetxt(path_bins + os.sep + names[i], (l_save,dl_save,dls_save))
	
	data = np.array([l,dl,dls]).T
	with open(name_file,'w+') as write_data:
		np.savetxt(write_data, data, fmt=['%f','%f','%f'])
	print "Tempo loop " +str(i+1)+":  " + str(time()-t)
np.savetxt(path_bins+os.sep+'z.dat',zz)
print "\nArquivos com redshifts criado."	
#load = np.loadtxt(path_bins + os.sep + names[0])
