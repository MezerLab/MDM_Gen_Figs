function fig_7_supfig_19_20(main_dir)

% this function generates figure 7: The “mosaic” nature of molecular and volumetric aging trajectories

%% Bring data 

% load human brain MRI data
load(fullfile(main_dir,'/fig7/human_data_fig7.mat'));
% this .mat file contain:
% - diff_AllMarakers is a 14X4 matrix with the aging effect-size in 14 ROIs
%   for 4 different aging markers {'Chemophysical','Water','Atrophy','Iron'}
% - std_AllMarakers is similar to diff_AllMarakers but with error
%   information.
% - diff_MDM is a 1X14X4 matrix with the aging effect size in 4 different
%   MDM dimentions {'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'} of 14
%   ROIs for 

%% generate figure 7

% figure 7A: Distinct spatial patterns of different aging markers throughout the brain
fig_7A(diff_AllMarkers,std_AllMarkers);

% figure 7B: The spatial correlation between aging-related changes in different biological markers. 
fig_7B(diff_AllMarkers);

%% generate supfigs 19-20

% sup. figure 19A: Distinct spatial patterns of different aging markers throughout the brain
supfig19A(diff_AllMarkers,diff_MDM);

% sup. figure 19B:  The spatial correlation between aging-related changes in volume, iron and water to all MDM dimensions. 
supfig19B(diff_MDM,diff_AllMarkers);

% sup. figure 20:  Distinct molecular aging trajectories- comparison of different MDM dimensions. 
supfig20(diff_MDM);

end