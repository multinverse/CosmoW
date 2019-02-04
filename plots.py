import numpy as np
import matplotlib.pyplot as plt
import Limber_approx as li
import parameters as par
from matplotlib import rc
from matplotlib.ticker import (NullFormatter,LogFormatter)  # useful for `logit` scale
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


majorLocator = MultipleLocator(20)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(5)

#from matplotlib.ticker import MultipleLocator
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
#rc('font', family='serif')

plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode'] = True



#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================

k = np.logspace(-4,np.log10(1), num=1000)
l = np.arange(1,1001,1)
z_Pk = 0
z_Cl = [0.127,0.138]

plot = {'pk':False,'cl':True,'pk_bao':False,'cl_bao':False}
params = {"R":par.R_sigma,"sigma_R":par.sigma8,"Om_h2":par.Om_h2,"Ob_h2":par.Ob_h2,"Oc_h2":par.Oc_h2,"Or_h2":par.Or_h2,
          "Og_h2":par.Og_h2, "Ok_h2":par.Ok_h2, "Od_h2":par.Od_h2,"h":par.h,"TCMB":par.TCMB,"ns":par.ns,'w':par.w,'model_w':par.model_w}
info_plot = {'k':k,'l':l,'eisenstein':True, 'nobao':True, 'nobaryon':True, 'bbks':True}
info_spec = {"spectrum": "eisenstein","amplitude": "cobe", "bias":"constant","OmegaHI_model":"constant","window":"battye"} #info_plot = 'BAO'#'single',smooth'


#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================



if plot['pk']==True:
	if info_plot['eisenstein'] == True:
		info_spec['spectrum'] = 'eisenstein'
		Pk = li.Pk_HI(k,z_Pk,params,info_spec)
		plt.loglog(k,Pk, label=r'$Eisenstein\&Hu$')
	if info_plot['nobao'] == True:
		info_spec['spectrum'] = 'nobao'
		Pk = li.Pk_HI(k,z_Pk,params,info_spec)
		plt.loglog(k,Pk, label = r'$No\, BAO$')		
	if info_plot['nobaryon'] == True:
		info_spec['spectrum'] = 'nobaryon'
		Pk = li.Pk_HI(k,z_Pk,params,info_spec)
		plt.loglog(k,Pk, label = r'$No\, Baryon$')
	if info_plot['bbks'] == True:
		info_spec['spectrum'] = 'bbks'
		Pk = li.Pk_HI(k,z_Pk,params,info_spec)
		plt.loglog(k,Pk, label = r'$BBKS$' )
	
	plt.xlim(np.amin(k),np.amax(k))
	leg = plt.legend(loc='best', shadow=False, fancybox=True)
	leg.get_frame().set_alpha(0)
	plt.xlabel(r'$k\,\left(Mpc^{-1}\right)$')
	plt.ylabel(r"$P\left(k\right)\,\left(Mpc\right)^3$")
	plt.show()



if plot['cl']==True:
	fig, ax = plt.subplots(2, 1) #(lins=..,columns=..,...)
	
	if info_plot['eisenstein'] == True:
		info_spec['spectrum'] = 'eisenstein'
		dl = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
		ax[0].loglog(l,dl, label=r'$Eisenstein\&Hu$')
	#	plt.xticks(np.arange(np.amin(l),np.amax(l),10))
	if info_plot['nobao'] == True:
		info_spec['spectrum'] = 'nobao'
		dls = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
		ax[0].loglog(l,dls,  label = r'$No\, BAO$')		
			

	ax[0].set_xlim(np.amin(l),np.amax(l))
	ax[0].set_xlabel(r'$\ell$')
	ax[0].set_ylabel(r"$D_{\ell}\,\left(\mu K\right)^3$")
	
	if info_plot['nobao'] == True and info_plot['eisenstein'] == True:
		ax[1].set_xlim(21,100)
		ax[1].plot(l,dl,  label = r'$Eisenstein\&Hu$')
		ax[1].plot(l,dls,  label = r'$No\, BAO$')
		ax[1].set_xlabel(r'$\ell$')
		ax[1].set_ylabel(r"$D_{\ell}\,\left(\mu K\right)^3$")
		ax[1].set_xticks(np.arange(20,101,10))
	plt.tight_layout()
	plt.show()
'''

if plot['cl']==True:
	if info_plot['eisenstein'] == True:
		info_spec['spectrum'] = 'eisenstein'
		dl = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
	if info_plot['nobao'] == True:
		info_spec['spectrum'] = 'nobao'
		dls = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
			
	if info_plot['nobao'] == True and info_plot['eisenstein'] == True:
		plt.xlim(np.amin(l),np.amax(l))
		#plt.plot(l,dl,  label = r'$Eisenstein\&Hu$')
		#plt.plot(l,dls,  label = r'$No\, BAO$')
		plt.loglog(l,dls,  label = r'$No\, BAO$')
		plt.loglog(l,dl,  label = r'$Eisenstein\&Hu$')
		plt.xlabel(r'$\ell$')
		plt.ylabel(r"$D_{\ell}\,\left(\mu K\right)^3$")
		plt.xticks(np.arange(21,100,10))
	plt.tight_layout()
	plt.show()
'''

if plot['pk_bao']==True:
		info_spec['spectrum'] = 'eisenstein'
		Pk = li.Pk_HI(k,z_Pk,params,info_spec)
		info_spec['spectrum'] = 'nobao'
		Pks = li.Pk_HI(k,z_Pk,params,info_spec)
		plt.subplot(211)
		plt.semilogx(k,Pk/Pks)
		#plt.xlim(1,200)
		plt.axhline(y=1.,color = "black")
		plt.xlim(np.amin(k),np.amax(k))
		plt.xlabel(r'$k\,\left(Mpc^{-1}\right)$')
		plt.ylabel(r"$P\left(k\right)/P_s\left(k\right)$")
		plt.subplot(212)
		plt.semilogx(k,Pk/Pks)
		plt.axhline(y=1.,color = "black")
		plt.xlabel(r'$k\,\left(Mpc^{-1}\right)$')
		plt.ylabel(r"$P\left(k\right)/P_s\left(k\right)$")
		plt.xlim(0.02,0.2)
		plt.show()
	
if plot['cl_bao']==True:
		info_spec['spectrum'] = 'eisenstein'
		dl = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
		info_spec['spectrum'] = 'nobao'
		dls = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
		plt.subplot(211)
		plt.semilogx(l,dl/dls)
		plt.xlim(1,200)
		plt.axhline(y=1.,color = "black")
		
		plt.ylabel(r"$C_{\ell}/C_{\ell,s}$")
		plt.subplot(212)
		plt.semilogx(l,dl/dls)
		plt.axhline(y=1.,color = "black")
		plt.xlabel(r'$\ell$')
		plt.ylabel(r"$C_{\ell}/C_{\ell,s}$")
		plt.xlim(21,100)
		plt.show()



if plot['cl_bao']==True:
		info_spec['spectrum'] = 'eisenstein'
		dl = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
		info_spec['spectrum'] = 'nobao'
		dls = li.Dl_HI_Limber(l,z_Cl,params, info_spec)
		err       = np.sqrt(2./(2.*l+1.))*(dl/dls)
		plt.errorbar(l, dl/dls, yerr=err, fmt='.')
		plt.axhline(y=1.,color = "black")
		plt.ylabel(r"$C_{\ell}/C_{\ell,s}$")
		plt.semilogx(l,dl/dls)
		plt.xlim(21,100)
		plt.show()
