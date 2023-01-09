#Root file with histograms from time-sliced fluxes for FD

ND_Laura_sliced_flux.root 



#Root file with histograms from time-sliced fluxes for FD

FD_Laura_sliced_flux.root


#Macro to plot fluxes

testflux.C


# nuTree2GlobesFlux.C & nuTree2GlobesFlux.h
nuTree2GlobesFlux convolves the relative arrival time of the neutrinos with the beam time structure and a detector resolution and makes histograms/globes flux files binned by arrival-time bins. The input to this code is : all_flux_laura_v3r5p4_nu_dec2422.root  
