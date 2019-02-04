Ob_h2         = 0.02230      
Oc_h2         = 0.1188       
Og_h2         = 2.472906693752628*1.e-5
Ok_h2         = 0.   
ns            = 0.9667       
As            = 2.143*1.e-9
Hubble_0      = 67.74        
TCMB          = 2.72548
sigma8        = 0.8159   #Planck2015
sigma8_camb   = 0.8071   #CAMB
Neff          = 3.046


h             = (Hubble_0/100.)
h2            = h**2
Om_h2         = Ob_h2+Oc_h2
R_sigma       = 8./h     #8 h-1Mpc

kmax          = 10   #Mpc-1
lmin          = 21
lmax          = 100 

def Omega_radh2(Ogammah2 = Og_h2,neff=Neff): return Ogammah2*(1. + 0.2271*neff)
Or_h2    = Omega_radh2()
Od_h2    = h2 - (Om_h2 + Ok_h2 + Or_h2)
model_w  = 'const'       # "const","cpl","EDE"
w        = -1.           # -1.019
wa       =  0            #if model = "cpl"


#z_coinc = -1. + ((Od_h2/h2/(1. - Od_h2/h2))**(1./3.)), -1. + ((Od_h2/Om_h2)**(1./3.)) 

#params    = {"R":R_sigma,"sigma_R":sigma8,"Om_h2":Om_h2,"Ob_h2":Ob_h2,"Oc_h2":Oc_h2,"Or_h2":Or_h2,"Og_h2":Og_h2,"Ok_h2":Ok_h2,"Od_h2":Od_h2,"h":h,"TCMB":TCMB,"ns":ns,"w":w,"model_w":model_w}

#Para plotar os espectros
#plot      = {"pk":False,"cl":False,"pk_bao":True,"cl_bao":False}
#info_plot = {'k':k,'l':l,'eisenstein':True, 'nobao':True, 'nobaryon':True, 'bbks':True}
#info_spec = {"spectrum": "eisenstein","amplitude": "cobe", "bias":"constant","OmegaHI_model":"constant","window":"battye"}


#Para fazer MCMC dos fits
