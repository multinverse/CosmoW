#! /usr/bin/env python
# -*- coding: utf-8 -*-

#import MCMC_parameters as mpar

########################################################################
#
# Define the configurations of the GetDist
#
########################################################################

mylabels   = [r"$A$", r"$\alpha_{\perp}$",r"$\Sigma$"]
def mylegends(zeff):return [str(int(10000*zeff)/10000.)]

plot_triangle          = True
plot_2d                = True
print_terminal_datas   = True
save_datas             = True
shaded_getdist_plot    = False
contour_filled         = True
getdist_sigma_error_bf = 1                              #1 (default) para regiao de 68%, limit=2 para a segunda regiao de contorno (especificada pelo  self.contours). VER: http://getdist.readthedocs.io/en/latest/mcsamples.html, GetInlineLatex().
		
