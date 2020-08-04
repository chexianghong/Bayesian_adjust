# Bayesian_adjust
the code of making Landsat consistent with MODIS NBAR using Bayesian theory 

Step1:
  run produce_sigma.py to get the optimal sigma for MODIS NBAR PSF 
  
Step2:
  run bayesian_adjust.py to adjust Landsat data to obtain fully MODIS-consistent Landsat data
  
Step3:
  run rma_adjust.py to obtain the linearly MODIS-consistent Landsat data
