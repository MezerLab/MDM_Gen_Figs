
function fig3_supfig9(main_dir,R2s_correction)

% this function generates Fig. 3: The biological interpretation of the MDM signatures based on comparison between in-vivo and post mortem data.

%% bring data

% load human MRI data
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

% if R2s_correction==1 use data that was corrected for R2*:
if R2s_correction==1
    load(fullfile(main_dir,'/supfig_21_29/human_brain_data_R2s_corrected.mat')); 
end

% bring histological data from "Lipid compositions of different regions of
% the human brain during aging" (Söderberg 1990):
lipids=gen_lipids_mat();
LipidsMat=lipids.percent;
 
set(0, 'DefaultFigureVisible', 'on');

%% arragne MRI data

Coef=[1,3,4]; % MDM / MTV / standard qMRI parameters 
MRIMat=nan(length(lipids.FS),4,3); % this is a matrix of 7 ROIs X qMRI parameters X MDM  / MTV / standard qMRI parameters 
MRIstd=nan(length(lipids.FS),4,3); % same as MRImat but for the errors


% MTV maps of four subjects had bias in the lower part of the brain and they
% were therefore excluded from the analysis presented in fig. 3 which includes
% ROIs in the brainstem:
young_ind([7,13,16,21])=[];

% Average MRI data and average over subjects:
allsub=nan(length(lipids.FS),46,length(fit_info.str),3);
for jj=1:4
    for k=1:length(Coef) % MDM / MTV / standard qMRI parameters 
        par=squeeze(fit_info.data{jj}(:,Coef(k),young_ind));
        for ii=1:length(lipids.FS) % loop over ROIs
            tmp=[];
            tmp=par(lipids.FS{ii},:);
            MRIMat(ii,jj,k)=nanmean(tmp(:));
            MRIstd(ii,jj,k)=nanstd(tmp(:));
            allsub(ii,1:length(tmp(:)),jj,k)=tmp(:);
        end
    end
end

%% color
color={[  18 93 119]./255,[34 100 95]./255,[1 153 71]./255,[135,190,36]./255,[  151 67 20]./255, [153 120 175]./255, [ 216 120 147]./255};
color=color(1:length(lipids.ROI));
color=cell2mat(color);
color=reshape(color,[3,length(lipids.ROI)]);
color=color';

%% Plot Fig. 3A-B: 

% compute the 1st PC of ex-vivo molecular compusition:
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(LipidsMat(:,[1,2,3,4,5,6,7])));
% project error on the PC:
red_err=lipids.err(:,[1,2,3,4,5,6,7])./std(LipidsMat(:,[1,2,3,4,5,6,7]));
red_err=red_err*coeff;

% compute the 1st PC of in-vivo MDM:
[coeffMRI,scoreMRI,latent,tsquared,explainedMRI,mu] = pca(zscore(MRIMat(:,[1,2,3,4],1)));

% for R2* correction project corrected data on the PC of the original
% figure:
if R2s_correction==1
coeff_new=[0.6680 ,     0.6701   ,0.0239, -0.3227 ]';
scoreMRI=zscore(MRIMat(:,[1,2,3,4],1))*coeff_new;
end

% project error on the PC:
err_MRI=MRIstd(:,[1,2,3,4],1)./std(MRIMat(:,[1,2,3,4],1));
err_MRI=err_MRI*coeffMRI;

% plot figure 3B: The similarity between the ex-vivo molecular variability and the in-vivo MDM variability across brain regions.
pc=1;
pcMRI=1;
h=[];
s=figure;
subplot(1,2,2)
hold on
% plot data points with errorbars:
for r=1:length(lipids.ROI)
    h(r)=scatter(scoreMRI(r,pcMRI),score(r,pc),100,color(r,:),'filled');
    errorbar(scoreMRI(r,pcMRI),score(r,pc),abs(red_err(r,pc)),'Color',color(r,:));
    errorbar(scoreMRI(r,pcMRI),score(r,pc),abs(err_MRI(r,pcMRI)),'horizontal','Color',color(r,:));
end
% evaluate the correlation:
mdl = fitlm(scoreMRI(:,pcMRI),score(:,pc));
linvec=[-3,3];
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
% plot fitted line:
plot(linvec,yprediction,'k--');
xlim([-3 3]);
grid on
hold off
xlabel('1st PC in-vivo MDM');
ylabel('1st PC ex-vivo Molecular');
legend(h,lipids.ROI,'Location','southeast','FontSize',13);
% plot significance:
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f')),[stars]};
xlim([-2 3]);
ylim([-3 3]);
set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18
annotation('textbox',[0.57267 0.76635 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14)

% plot figure 3A: Establishing the agreement between the postmortem dataset and the in vivo MRI measurements.
subplot(1,2,1);
hold on
% plot data points with errorbars:
for r=1:length(lipids.ROI)
    h(r)=scatter(MRIMat(r,2,2),lipids.TotalPhos(r),100,color(r,:),'filled');
    errorbar(MRIMat(r,2,2),lipids.TotalPhos(r),MRIstd(r,1,2),'horizontal','Color',color(r,:));
    errorbar(MRIMat(r,2,2),lipids.TotalPhos(r),lipids.TotalPhoserr(r),'vertical','Color',color(r,:));
end
% evaluate the correlation:
mdl = fitlm(MRIMat(:,2,2),lipids.TotalPhos);
linvec=[0.1 0.39];
xlim([0.15 0.39]);
ylim([15 40]);
set(gca, 'ytick',[15:5:40]);
% plot fitted line:
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
plot(linvec,yprediction,'k--');
grid on
hold off
xlabel('MTV [fraction]');
ylabel('Total Phospholipid [mg/g wet weight]');
% plot significance:
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f')),[ stars]};
annotation('textbox',[0.13255 0.76635 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14)
set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18
set(gcf, 'Position',[1 1 1268 500]);


%% Fig. 3C: Predicting molecular composition with MRI. 

lab={'PE/PdtCho [fraction]','Spg [p.u.]','Phospholipids/Proteins [fraction]'};
laby={'Prediction [fraction]','Prediction [p.u.]','Prediction [fraction]'};

% try to predict lipids that account for most of the molecular variability:
LipidsMatnew=[];
LipidsMatnew(:,1)=LipidsMat(:,1)./LipidsMat(:,3); % compute the PE/PtdCho ratio
LipidsMatnew(:,[2,3])=LipidsMat(:,[5 7]); % Spg and the ratio of phopholipids/proteins

R=[];
P=[];
prediction=[];

% compute cross validation prediction from MDM to lipid composition:

A=(MRIMat(:,[ 1 2 ],1)); % MDM measurements for dR1/dMTV and dMTsat/dMTV
B=LipidsMatnew(:,:); % lipid composition

% solve linear equation:
for rr=1:length(lipids.ROI)
    X = linsolve(A([1:rr-1,rr+1:end],:),B([1:rr-1,rr+1:end],:));
    prediction(rr,:)=A(rr,:)*X;
end

% evaluate the prediction:
for ii=1:3
    mdl = fitlm(B(:,ii),prediction(:,ii));
    R(ii)=mdl.Rsquared.Adjusted;
    P(ii)=mdl.Coefficients{2,4};
end

% plot figure:
s=figure('Position', [583 489 1110 450]);
for ii=1:size(LipidsMatnew,2)
    subplot(1,size(LipidsMatnew,2),ii)
    hold on
    % plot data points:
    for r=1:length(lipids.ROI)
        h(r)=scatter(LipidsMatnew(r,ii),prediction(r,ii),200,color(r,:),'filled');
    end
    xlabel([lab{ii}]);
    ylabel([laby{ii}]);
    hold off
    grid on
    if ii==1
        xlim([0.5 2.2]);
        ylim([0.5 2.2]);
    elseif ii==2
        xlim([5.7 19.5]);
        ylim([5.7 19.5]);
        set(gca, 'xtick',[6:4:18]);
        
    elseif ii==3
        xlim([0.14 0.42]);
        ylim([0.14 0.42]);
    end
    % plot identity line:
    hline = refline(1,0);
    hline.Color = 'k';
    hline.LineStyle='--';
    hold off
    set(gca,'FontSize',15);
    % plot significance:
    [stars]=pval2stars(P(ii),'num');
    eq={strcat('R^2=',num2str(R(ii),'%0.2f')),[stars]};
    ypos=max(get(gca, 'ylim'));
    xpos=min(get(gca, 'xlim'));
    text(xpos+0.1*xpos,ypos-0.05*ypos,eq,'FontSize',15)
    
end
set(gcf,'color','w');
set(gcf, 'Position',[1 1 1268*1.5 500]);


%% SupFig. 9A: The correlation of individual lipids (derived from the literature) with the 1st principal component (PC) of in-vivo MDM signatures. 

% compute the 1st PC of in-vivo MDM:
[coeffMRI,scoreMRI,latent,tsquared,explainedMRI,mu] = pca(zscore(MRIMat(:,[1,2,3,4],1)));
% project errors on PC:
err_MRI=MRIstd(:,[1,2,3,4],1)./std(MRIMat(:,[1,2,3,4],1));
err_MRI=err_MRI*coeffMRI;

% plot figure
pcMRI=1;
h=[];
s=figure;
for pc=1:length(lipids.all)
    subplot(3,3,pc)
    hold on
    for r=1:length(lipids.ROI)
        h(r)=scatter(scoreMRI(r,pcMRI),LipidsMat(r,pc),100,color(r,:),'filled');
        errorbar(scoreMRI(r,pcMRI),LipidsMat(r,pc),abs(lipids.err(r,pc)),'Color',color(r,:));
        errorbar(scoreMRI(r,pcMRI),LipidsMat(r,pc),abs(err_MRI(r,pcMRI)),'horizontal','Color',color(r,:));
    end
    % evaluate the correlation:
    mdl = fitlm(scoreMRI(:,pcMRI),LipidsMat(:,pc));
    linvec=[-3,3];
    yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
    % plot fitted line:
    plot(linvec,yprediction,'k--');
    xlim([-3 3]);
    grid on
    hold off
    xlabel('1st PC in-vivo MDM');
    ylabel(lipids.all{pc});
    % plot significance:
    [stars]=pval2stars(mdl.Coefficients{2,4},'num');
    eq=['R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f'),' ',[stars]];
    title(eq,'FontWeight','Normal');
    xlim([-2 3]);
    set(gca,'FontSize',15) 
end
legend(h,lipids.ROI,'Location','southeast','FontSize',13);
set(gcf, 'Position',[1 1 1250 1050]);

%% SupFig. 9B: The correlation of the molecular variability with standard qMRI parameters is lower compared to the correlation with MDM (Fig. 3B).

% compute the 1st PC of ex-vivo molecular composition:
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(LipidsMat(:,[1,2,3,4,5,6,7])));
% project errors on PC:
red_err=lipids.err(:,[1,2,3,4,5,6,7])./std(LipidsMat(:,[1,2,3,4,5,6,7]));
red_err=red_err*coeff;

% compute the 1st PC of in-vivo standard qMRI parameters:
[coeffMRI,scoreMRI,latent,tsquared,explainedMRI,mu] = pca(zscore(MRIMat(:,[1,2,3,4],3)));
% project errors on PC:
err_MRI=MRIstd(:,[1,2,3,4],3)./std(MRIMat(:,[1,2,3,4],3));
err_MRI=err_MRI*coeffMRI;

% plot figure:
pc=1;
pcMRI=1;
h=[];
s=figure;
subplot(1,2,1)
hold on
% plot data points and errorbars:
for r=1:length(lipids.ROI) 
    h(r)=scatter(scoreMRI(r,pcMRI),score(r,pc),100,color(r,:),'filled');
    errorbar(scoreMRI(r,pcMRI),score(r,pc),abs(red_err(r,pc)),'Color',color(r,:));
    errorbar(scoreMRI(r,pcMRI),score(r,pc),abs(err_MRI(r,pcMRI)),'horizontal','Color',color(r,:));
    
end
% evaluate the correlation:
mdl = fitlm(scoreMRI(:,pcMRI),score(:,pc));
linvec=[-3,3];
xlim([-3 3]);
% plot fitted line:
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
plot(linvec,yprediction,'k--');
grid on
hold off
xlabel('1st PC in-vivo MRI standard');
ylabel('1st PC ex-vivo Molecular');
% plot significance:
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f')),[ stars]};
set(gca,'FontSize',15) 
annotation('textbox',[0.13255 0.76635 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14)

%% SupFig. 9C: The correlation of the molecular variability with MTV is lower compared to the correlation with MDM (Fig. 3B)

% compute the 1st PC of ex-vivo molecular composition:
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(LipidsMat(:,[1,2,3,4,5,6,7])));
% project errors on PC:
red_err=lipids.err(:,[1,2,3,4,5,6,7])./std(LipidsMat(:,[1,2,3,4,5,6,7]));
red_err=red_err*coeff;

% compare to MTV:
subplot(1,2,2)
hold on

% plot data points and errorbars:
for r=1:length(lipids.ROI)
    h(r)=scatter(MRIMat(r,2,2),score(r,pc),100,color(r,:),'filled');
    errorbar(MRIMat(r,2,2),score(r,pc),MRIstd(r,1,2),'horizontal','Color',color(r,:));
    errorbar(MRIMat(r,2,2),score(r,pc),abs(red_err(r,pc)),'vertical','Color',color(r,:));
end
% evaluate the correlation:
mdl = fitlm(MRIMat(:,2,2),score(:,pc));
linvec=[0.1 0.39];
xlim([0.1 0.39]);
% plot fitted line:
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
plot(linvec,yprediction,'k--');
grid on
hold off
xlabel('MTV [fraction]');
ylabel('1st PC ex-vivo Molecular');
legend(h,lipids.ROI,'Location','southeast');
% plot significance:
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f')),[ stars]};
set(gca,'FontSize',15) 
set(gcf, 'Position',[1 1 1268 500]);
annotation('textbox',[0.57267 0.76635 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14)

end
