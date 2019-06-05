function fig1(main_dir)

% this function generates figure 1: The dependency of qMRI parameters on the molecular composition in lipids samples. 

% Fig. 1A-C and SupFig. 1A-B: The dependency of R1, R2 and MTsat on MTV for different lipids.  
fig1ABC(main_dir)

% Fig 1D and SupFig. 1C: Predicting the MRI signal of a lipid mixture from the signal of pure lipids. 
fig1D(main_dir)

% Fig 1E: Predicting the lipid composition of 12 mixtures using the MDM method.
fig1E(main_dir)

% SupFig. 1D: The dependency of R1, R2 and MTsat on MTV for non-lipid compound  
supfig1D(main_dir)

% SupFig. 2: MDM-based prediction for lipids samples with cross validation. 
supfig2(main_dir) 

end
