function fig2_supfig_3_8(main_dir)

% this fucntion generates Figure 2: The MDM method provides region-specific signatures in the in-vivo human brain. 
% and sup. figures 3-8
%% bring data

% load human data
load(fullfile(main_dir,'/fig2_5_6/human_data_1.mat'));
% this .mat file contains:
% - fit_info is a the structure of MDM measurements for all subjects. For each qMRI parameter (fit_info.str) the
%   measurements are stored in the data field in a 35X6X41 matrices. There
%   are 35 ROIs and 41 subjects. The 2nd dimention represent: [slope,
%   intersection, mean MTV, mean qMRI parameter, R^2 and cvRMSE] for each
%   ROI.
% - old_ind young_ind indicates the age group of each subject.
% - ROI_list represents the names of the different ROIs.
% - fit_infoCTX is similar to fit_info but for cortical ROIs
% - C represents the freesurfer labels of the cortical ROIs.

%% settings

v=[1:5,8:16]; % left label of ROIs for the analysis
r=[20:24,27:35]; % right label of ROIs for the analysis
l=v; % left label of ROIs for the analysis

comrl=1; % join left and right ROIs if 1.

if comrl % fix labels when right and left ROIs are joined
    ROI_list(l)=cellfun(@(x) x(6:end),ROI_list(l),'UniformOutput',false);
end

% New Colormap
color_new={[106 54 136 ]./255; ...
    [153 120 175]./255;...
    [203 177 218]./255;...
    [ 214 86 121]./255;...
    [ 216 120 147]./255;...
    [  244 173 199]./255;...
    [  151 67 20]./255;...
    [  215 102 0]./255;...
    [  246 188 65]./255;...
    [  218 127 2]./255;...
    [  18 93 119]./255;...
    [ 0 133 163]./255;...
    [ 133 198 220]./255;...
    [ 5 149 198]./255};

set(0, 'DefaultFigureVisible', 'on');
%% figure 2 BC: Calculating the MDM signatures. 

fig2BC(main_dir,ROI_list,fit_info,v,l,r,color_new,comrl,young_ind)


%% figure 2D: Unique MDM signatures for different brain regions

fig2d_spider(color_new,comrl,young_ind,fit_info,v,l,r)


%% plot sup. figures 3-8:
supfig_3_8(main_dir,ROI_list,fit_info,v,l,r,color_new,comrl,young_ind)

end