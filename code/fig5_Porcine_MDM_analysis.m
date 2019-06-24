function fig5_Porcine_MDM_analysis(main_dir)

% this function generates figure 5: Post mortem validation of the MDM approach

%% Load MRI and histology data of two porcine brains

data=load(fullfile(main_dir,'/fig4/porcine_brain_data.mat'));
ROIs_labels=data.ROIs_labels; % ROIs names
MDM=data.MDM; % ROIsX [slope,intersection, mean MTV, mean qMRI parameter, R^2 and cvRMSE]X qMRI parameters (MTsat,R1,MD.R2)
MTV_MRI=data.MTV_MRI; % MTV values from MRI for each ROI
MTV_histology=data.MTV_histology;  % MTV values from histology for each ROI
lipids=data.lipids; % lipid compositon (7 lipids)  for 30 ROIs
lipids_labels=data.lipids_labels; % names of the detected lipids

%% Plot figure 5: Post mortem validation of the MDM approach

color=[0.03,0.04,0.35];

% plot figure 5F: The dependency of MTsat on MTV in three example brain regions. 
figure;
subplot(3,2,2)
fig4F(main_dir)

% plot figure 5G: Unique MDM signatures for different brain regions. 
subplot(3,2,4)
fig4G(MDM);

% plot figure 5H: The MDM signatures explain the molecular variability across brain regions. 
subplot(3,2,6)
[mdl,slopes_score1,score1]=compare_lipids_post_mortem(ROIs_labels,MDM,lipids./10^9,lipids_labels,1,0,color);

% plot figure 5C: Validation of the qMRI estimation of MTV.  
subplot(3,2,1)
[mdl]=compare_mtv_post_mortem_allbrains(MTV_histology,MTV_MRI,color);

% plot figure 5E: Standard qMRI parameters do not explain the molecular variability
subplot(3,2,5)
[mdl,slopes_score4,score4]=compare_lipids_post_mortem(ROIs_labels,MDM,lipids./10^9,lipids_labels,4,0,color);

% plot figure 5D: MTV does not explain the molecular variability
subplot(3,2,3)
[mdl,slopes_score3,score3]=compare_lipids_post_mortem(ROIs_labels,MDM,lipids./10^9,lipids_labels,3,0,color);

set(gcf,'color','white')
set(gcf, 'Position',[1 1 1100 1000*1.5]);

%% Sup. Fig. 10: Histological and MDM analysis of the porcine brain and comparison to human. 

% plot Sup. Figure 10A:  Phospholipid composition of the porcine and human brain. 

% bring human brain lipid compositon (from the analysis of figure 3):
lipids2=gen_lipids_mat();
LipidsMat=lipids2.percent;
LipidsName=lipids2.all;

polarlipid=lipids(:,[3,4,6,7]); % show only polar lipids (as they are shown in fig. 3)
figure
subplot(1,2,1)
% pie chart of the median of each lipid across all porcine ROIs 
p=pie(nanmedian(polarlipid,1));
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = LipidsName([3,1,2,5]);
pText(1).String = [txt{1}(1:end-7),': ',percentValues{1}];
pText(2).String = [txt{2}(1:end-7),': ',percentValues{2}];
pText(3).String = [txt{3}(1:end-7),': ',percentValues{3}];
pText(4).String = [txt{4}(1:end-7),': ',percentValues{4}];
set(gca,'FontSize',15);

subplot(1,2,2)
% pie chart of the median of each lipid across all human ROIs 
p=pie(nanmedian(LipidsMat(:,[3,1,2,5]),1));
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = LipidsName([3,1,2,5]);
pText(1).String =[txt{1}(1:end-7),': ',percentValues{1}];
pText(2).String = [txt{2}(1:end-7),': ',percentValues{2}];
pText(3).String = [txt{3}(1:end-7),': ',percentValues{3}];
pText(4).String = [txt{4}(1:end-7),': ',percentValues{4}];
set(gca,'FontSize',15);

roi=1:30;
color_new=parula(30);

% plot sup. fig. 10B: Unique MDM signatures of 30 brain regions in two porcine brains. 
figure;
spiderpostmortem(MDM,roi,color_new)

% plot sup. fig. 10C: The similarity between the molecular variability and the MDM variability in the same brain after removing outliers. 
figure;
subplot(1,3,1)
[mdl]=compare_outliers_post_mortem(color,slopes_score1,score1,1);
subplot(1,3,2)
[mdl]=compare_outliers_post_mortem(color,slopes_score4,score4,4);
subplot(1,3,3)
[mdl]=compare_outliers_post_mortem(color,slopes_score3,score3,3);

set(gcf, 'Position',[1 1 1500 1000*0.5]);
end