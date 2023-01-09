#Root file with histograms from time-sliced fluxes for FD

ND_Laura_sliced_flux.root 



#Root file with histograms from time-sliced fluxes for FD

FD_Laura_sliced_flux.root


#Macro to plot fluxes

testflux.C


# nuTree2GlobesFlux.C & nuTree2GlobesFlux.h
nuTree2GlobesFlux convolves the relative arrival time of the neutrinos with the beam time structure and a detector resolution and makes histograms/globes flux files binned by arrival-time bins. The inputs to this code are : all_flux_laura_v3r5p4_nu_dec2422.root, 531mhz_down10_up20_buff30.root. 

#To run the code, in a ROOT session, you can do:

Root > .L nuTree2GlobesFlux.C

Root > nuTree2GlobesFlux t

Root > t.Loop();   // Loop on all entries

The input file all_flux_laura_v3r5p4_nu_dec2422.root can be downloded from here

https://drive.google.com/file/d/1KN5EyuSzy9jctRaXLBXkoGWgEfY94BnQ/view?usp=share_link
