#! /usr/bin/env python
# -*- coding: utf-8 -*-

from math import pi
from parameters import h


#------------------------------
#21cm BINGO--------------------
#------------------------------
Lambda_21cm     = 0.21106     #m
freq_21cm       = 1420.405751 #MHz


n0              = 0.03*(h**3)  #sources comovel numerical density. [n0.h3] = Mpc-3
Omega_Survey    = 2900.    #deg2
V_survey        = 1.2      #Gpc3
Theta_FWHM      = 40.      #arcmin
V_pix           = 810.     #Mpc3
T_sys           = 50.      #Kelvin
n_f             = 50.      #Horns
t_obs           = 1.       #Observation time in years
z_max_BINGO     = 0.48
z_min_BINGO     = 0.13
freq_max_BINGO  = 1260.     #MHz
freq_min_BINGO  = 960.      #MHz
bin_freq_BINGO  = 7.5       #MHz  *Warning! This is to 40bins
D_antenna_diam  = 100.      #m,    see Battye et al. 2012, pag 9.


n_bins          = 15 #15,70,300 to have del_z=0.02,0.005,0.001   
del_freq        = (freq_max_BINGO-freq_min_BINGO)/n_bins         #MHz. bandwidth

gamma_max_BINGO = 2*15*pi/360.         #L.Olivari et. al. (arxiv: 1707.07647v3, pg:7,seção:4.2).
bin_l           = pi/gamma_max_BINGO   #Utilizar tal bin de multipolos para analise faz com que podemos desprezar as correlações entre os multipolos, e os bins são aproximadamente independentes.
