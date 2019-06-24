
% MDM_Gen_Figs
% ============
%
% Code associated with:
% Filo et. al. (2019). Disentangling molecular alterations from water-content changes in the aging
% human brain using quantitative MRI.
%     
% 
% This repository includes the code and data required to reporoduce all the
% figures in the paper. 
% To use this code:
% 1) Save this repository on your computer.
% 2) In the script MDM_Gen_Figs.m: update the local path of this repository in the variable 'main_dir'
% 3) Run the script MDM_Gen_Figs.m
% 
%
%  ## SOFTWARE REQUIREMENTS ###
%  * Code was tested on MATLAB R2017b
%  * boundedline-pkg - https://github.com/kakearney/boundedline-pkg 
%  * IoSR Matlab Toolbox - https://github.com/IoSR-Surrey/MatlabToolbox
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2019

%% Local path of the MDM_Gen_Figs repository

main_dir='/ems/elsc-labs/mezer-a/Mezer-Lab/projects/code/MTV_Vs_qMRI/MDM_Gen_Figs/';

%% Figure 1 and SupFig. 1-2:

fig1(main_dir);

%% Figure 2 and SupFig. 3-8:

fig2_supfig_3_8(main_dir)

%% Figure 3 and SupFig. 9

% figure 3 and sup. fig. 9: The biological interpretation of the MDM signatures based on comparison between in-vivo and post mortem data.
fig3_supfig9(main_dir,0)

%% Figure 4 and SupFig. 11-12

% figure 4 and sup. fig. 11-12: MDM correlation with specific gene expression patterns throughout the cortex. 
fig4_supfig11_12(main_dir)

%% Figure 5 and SupFig. 10: Post mortem validation of the MDM approach

fig5_Porcine_MDM_analysis(main_dir)

%% Figure 6-7 and SupFig. 13-18: Region-specific aging-related molecular changes revealed by the R1/MTsat dependency on MTV. 

fig_6_7_andSup_aging(main_dir)

%% Figure 8 and SupFig 19-20: The “mosaic” nature of molecular and volumetric aging trajectories

fig_8_supfig_19_20(main_dir)

%% Sup. Fig. 21-29:  R2* correction

supfig_21_29(main_dir)

%% Sup. Fig. 30:  Examination of the linear relationship of qMRI parameters and 1/PD. 
supfig30(main_dir)
