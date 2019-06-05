
function fig_5_6_andSup_aging(main_dir)

% This function generates Fig. 5-6: Region-specific aging-related molecular changes revealed by the R1/MTsat dependency on MTV. 
% and sup fig 13-18

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

% aarange ROI list:
ROI_list=strrep(ROI_list,'ctx','CTX');
ROI_list=strrep(ROI_list,'wm','WM');

str_vec=[1 2 3 4]; % which qMRI parameters to take
if length(fit_info.str)>4
    str_vec=[1 2 3 4 5];
end

%% seperate young and old data

% preforms a two-sample t-test to evaluate the significance of
% aging-related changes:
[hip,FDRval]=gen_young_old_matrix(young_ind,old_ind,fit_info,v,l,r,comrl,str_vec);
% hip is a matrix of [ROIs X age group X [slope,intersection,mean MTV, mean
% qMRI parameter] X qMRI parameters]
% FDRval is a matrix of FDR corrected p-values with dimentions of [slope,
% mean MTV, mean qMRI parameter] X ROIs X qMRI parameters.


% arrange data of young subjects for spider plots:
Test_young=nan(4,length(v),length(young_ind).*2); % qMRI parameters X ROIs X subjects
Test_young(:,:,1:length(young_ind))=[rot90(fit_info.data{1}(v,[1],young_ind)) ; rot90(fit_info.data{2}(v,[1],young_ind)); rot90(fit_info.data{3}(v,[1],young_ind)); rot90(fit_info.data{4}(v,[1],young_ind))];
for ii=1:length(v)
    if ismember(v(ii),l) && comrl
        ind=find(l==v(ii));
        Test_young(:,ii,length(young_ind)+1:end)=[rot90(fit_info.data{1}(r(ind),[1],young_ind)) ; rot90(fit_info.data{2}(r(ind),[1],young_ind)); rot90(fit_info.data{3}(r(ind),[1],young_ind)); rot90(fit_info.data{4}(r(ind),[1],young_ind))];
    end
end

% arrange data of old subjects for spider plots:
Test_old=nan(4,length(v),length(young_ind)*2); % qMRI parameters X ROIs X subjects
Test_old(:,:,1:length(old_ind))=[rot90(fit_info.data{1}(v,[1],old_ind)) ; rot90(fit_info.data{2}(v,[1],old_ind)); rot90(fit_info.data{3}(v,[1],old_ind)); rot90(fit_info.data{4}(v,[1],old_ind))];
for ii=1:length(v)
    if ismember(v(ii),l) && comrl
        ind=find(l==v(ii));
        Test_old(:,ii,length(old_ind)+1:length(old_ind)*2)=[rot90(fit_info.data{1}(r(ind),[1],old_ind)) ; rot90(fit_info.data{2}(r(ind),[1],old_ind)); rot90(fit_info.data{3}(r(ind),[1],old_ind)); rot90(fit_info.data{4}(r(ind),[1],old_ind))];
    end
end

color_spider=color_new;
color_spider=cell2mat(color_spider);

%% Fig. 5-6:  Region-specific aging-related molecular changes revealed by the R1/MTsat dependency on MTV. 

fig_5_6(hip,color_new,Test_young,Test_old,FDRval);
    
%% SupFig. 13: aging-related changes in MDM signatures

supfig13_spiderROIs(color_spider,FDRval,fit_info,Test_young,Test_old,v)

%% SupFig. 14-17: separating molecular and water related contributions. 

box_allROIs_fig(str_vec,hip,fit_info,ROI_list,color_new,FDRval)

%% SupFig. 18:   aging related changes revealed by the 1st principle component of MDM. 

box_allROIs_PC(hip,color_new,FDRval,ROI_list)

end

