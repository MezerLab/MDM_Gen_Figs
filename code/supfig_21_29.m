function supfig_21_29(main_dir)
% this function generates sup. fig. 21-29

% load human brain MRI data with R2* correction:
load(fullfile(main_dir,'/supfig_21_29/human_brain_data_R2s_corrected.mat'));

v=[1:5,8:16]; % ROIs for the analysis
r=[20:24,27:35]; % right label of ROIs for the analysis
l=v; % left label of ROIs for the analysis

str_vec=[1 2 3 4]; % which qMRI parameters to take
if length(fit_info.str)>4
    str_vec=[1 2 3 4 5];
end

comrl=1; % join left and right ROIs if 1.

if comrl % fix labels when right and left ROIs are joined
    ROI_list(l)=cellfun(@(x) x(6:end),ROI_list(l),'UniformOutput',false);
end

% arrange ROI list
ROI_list=strrep(ROI_list,'ctx','CTX');
ROI_list=strrep(ROI_list,'wm','WM');


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

%% Group yound and older subjects

% preforms a two-sample t-test to evaluate the significance of
% aging-related changes:
[hip,FDRval,R2mat]=gen_young_old_matrix(young_ind,old_ind,fit_info,v,l,r,comrl,str_vec);
% hip is a matrix of [ROIs X age group X [slope,intersection,mean MTV, mean
% qMRI parameter] X qMRI parameters]
% FDRval is a matrix of FDR corrected p-values with dimentions of [slope,
% mean MTV, mean qMRI parameter] X ROIs X qMRI parameters.
% R2mat is a matrix of R2* values for all ROIs- it is created only if there
% is R2* data in the fit_info matrix

%% sup. fig. 21: Iron-related changes between the two age groups. 

supfig21(R2mat,color_new,ROI_list,v,FDRval);

%% sup. fig. 22: Distinct regional dependencies on MTV following R2* correction. 

supfig22(ROI_list,v,fit_info,young_ind,l,comrl,r,color_new);

%%  supfig 23-24:  Aging-related changes revealed by the R1 dependency on MTV– R2* correction.

% Arragne data for spider plots:

% young subjects:
Test_young=nan(4,length(v),length(young_ind).*2); % qMRI parameter X ROIs X subjects
Test_young(:,:,1:length(young_ind))=[rot90(fit_info.data{1}(v,[1],young_ind)) ; rot90(fit_info.data{2}(v,[1],young_ind)); rot90(fit_info.data{3}(v,[1],young_ind)); rot90(fit_info.data{4}(v,[1],young_ind))];
for ii=1:length(v)
    if ismember(v(ii),l) && comrl
        ind=find(l==v(ii));
        Test_young(:,ii,length(young_ind)+1:end)=[rot90(fit_info.data{1}(r(ind),[1],young_ind)) ; rot90(fit_info.data{2}(r(ind),[1],young_ind)); rot90(fit_info.data{3}(r(ind),[1],young_ind)); rot90(fit_info.data{4}(r(ind),[1],young_ind))];
    end
end

% old subjects:
Test_old=nan(4,length(v),length(young_ind)*2); % qMRI parameter X ROIs X subjects
Test_old(:,:,1:length(old_ind))=[rot90(fit_info.data{1}(v,[1],old_ind)) ; rot90(fit_info.data{2}(v,[1],old_ind)); rot90(fit_info.data{3}(v,[1],old_ind)); rot90(fit_info.data{4}(v,[1],old_ind))];
for ii=1:length(v)
    if ismember(v(ii),l) && comrl
        ind=find(l==v(ii));
        Test_old(:,ii,length(old_ind)+1:length(old_ind)*2)=[rot90(fit_info.data{1}(r(ind),[1],old_ind)) ; rot90(fit_info.data{2}(r(ind),[1],old_ind)); rot90(fit_info.data{3}(r(ind),[1],old_ind)); rot90(fit_info.data{4}(r(ind),[1],old_ind))];
    end
end
color_spider=color_new;
color_spider=cell2mat(color_spider);

% generate figures:
fig_6_7(hip,color_new,Test_young,Test_old,FDRval);

%%  supfig 25-28: Aging-related changes in qMRI parameters revealed by the MDM approach – R2* correction. 

str_vec=1:4;
box_allROIs_fig(str_vec,hip,fit_info,ROI_list,color_new,FDRval)

%%  sup. fig. 29:  The biological interpretation of the MDM approach – following R2* correction. 

% figure 3:
fig3_supfig9(main_dir,1)

end