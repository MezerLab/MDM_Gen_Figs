function supfig11_12(main_dir)

% this function generates sup. fig. 11-12: MDM correlation with specific gene expression patterns throughout the cortex. 

%% load data

% load human MRI data
load(fullfile(main_dir,'/fig2_5_6/human_data_1.mat'));

% load gene expression labels
load(fullfile(main_dir,'/supfig11/genes_data.mat'));

% load gene expression data acquired from the Allen Brain Atlas
% (http://human.brain-map.org/well_data_files) and Ben-David & Shifman (2012):
load(fullfile(main_dir,'/supfig11/genes_data_2.mat'));
% each Seg_Data matrix has information of 84 ROIs. For each ROI you can
% find the eigengene of each of the 19 modules.
%% Arrange MRI data

subcorticalROI=[1:5,8,20:24,27]; % subcortical ROIs from freesurfer that are also in the gene expression data

corticalROI=1:68; % cortical ROIs in freesurfer

MRIpar=[1 2]; % work with the MTV derivatives of R1 and MTsat

% save MRI information in the following matrices:  ROIs X qMRI parameters X young subjects
total_slope_info=nan(length([subcorticalROI corticalROI]),length(MRIpar),length(young_ind)); 
total_normal_info=nan(length([subcorticalROI corticalROI]),length(MRIpar),length(young_ind));

% bring MDM measurements:
for ii=1:length(MRIpar)
    sub=fit_info.data{MRIpar(ii)}(subcorticalROI,1,young_ind);
    cortex=fit_infoCTX.data{MRIpar(ii)}(corticalROI,1,young_ind);
    total_slope_info(:,ii,:)=[sub; cortex];
    sub=fit_info.data{MRIpar(ii)}(subcorticalROI,4,young_ind);
    cortex=fit_infoCTX.data{MRIpar(ii)}(corticalROI,4,young_ind);
    total_normal_info(:,ii,:)=[sub; cortex];
end
total_MTV_info=[fit_info.data{2}(subcorticalROI,3,young_ind) ; fit_infoCTX.data{2}(corticalROI,3,young_ind)];

% average over subjects:
slope_median_data=nanmedian(total_slope_info,3);
normal_median_data=nanmedian(total_normal_info,3);
mtv_median_data=nanmedian(total_MTV_info,3);

% romove ROIs with no information:
exclude_regions=isnan(mean(slope_median_data,2));
slope_median_data=slope_median_data(~exclude_regions,:);
normal_median_data=normal_median_data(~exclude_regions,:);
mtv_median_data=mtv_median_data(~exclude_regions,:);

% arragne ROI list:
total_ROI = [fit_info.ROI(subcorticalROI); fit_infoCTX.ROI(corticalROI)];
total_ROI(exclude_regions)=[];
total_ROI_num=nan(length(total_ROI),1);
for ii=1:length(total_ROI)
    Index = find(strcmp(FS_label, total_ROI{ii}));
    total_ROI_num(ii)=FS_num(Index);
end

%% Arrange Gene Expression data

Median_Data=nan(length(total_ROI_num),19,2); % ROIs X modules X donors
for ii=1:length(total_ROI_num) % loop over ROIs
    Index = find(Seg_Data9861{:,1}==total_ROI_num(ii)); % find ROI in the gene expression table
    if ~isempty(Index)
        Median_Data(ii,:,1)=Seg_Data9861{Index,4}; % take the median eigengene of all modules for this ROI
        Index = find(Seg_Data10021{:,1}==total_ROI_num(ii)); 
        Median_Data(ii,:,2)=Seg_Data10021{Index,4}; % the same for the other donor
    else
        % if there is no data for this ROI place nans:
        Median_Data(ii,:,:)=nan(19,2);
        slope_median_data(ii,:)=nan(1,length(MRIpar));
        normal_median_data(ii,:)=nan(1,length(MRIpar));
        
    end
end
% avarage over the two donors:
Median_Data=nanmedian(Median_Data,3);

% remove ROIs with no gene expression data from the MRI information:
slope_median_data=slope_median_data(~isnan(Median_Data(:,1)),:);
normal_median_data=normal_median_data(~isnan(Median_Data(:,1)),:);
mtv_median_data=mtv_median_data(~isnan(Median_Data(:,1)),:);

% remove ROIs with no gene expression data from the gene expression information:
total_ROI=total_ROI(~isnan(Median_Data(:,1)),:);
total_ROI_num=total_ROI_num(~isnan(Median_Data(:,1)),:);
Median_Data=Median_Data(~isnan(Median_Data(:,1)),:);

%% Modules from Ben-David & Shifman (2012):

modules = {'Membrane','Synapse II','Mitochondrion','Synapse','Translational elongation','MEgreenyellow','MEblue','MEmagenta','MEcyan','Zinc ion binding','MElightyellow','MEblack','MEturquoise','MEgreen','Poly-Gly','MElightcyan','MEbrown','MEpink','MEpurple'};

%% Compare MRI and gene expression:

ROIvaec=[9:length(total_ROI)]; % take only cortical ROIs
par=1;

Gene=Median_Data(ROIvaec,:); % gene expression data

slopePCA=slope_median_data(ROIvaec,:); % MDM data
normalPCA=normal_median_data(ROIvaec,:); % standard qMRI data
mtv=mtv_median_data(ROIvaec,1); % MTV data

% compute PCA of MDM:
[slope_wcoeff,slope_score,slope_latent,slope_tsquared,slope_explained] = pca(zscore(slopePCA));

% compute PCA of standard qMRI parameters:
[normal_wcoeff,normal_score,normal_latent,normal_tsquared,normal_explained] = pca(zscore(normalPCA));

% join 1st PC of MDM, 1st PC of standard qMRI parameters and MTV:
MRIPar=[slope_score(:,1) normal_score(:,1) mtv];

% colors:
ModuleColor = 1/255*[169,169,169;144,238,144;255,0,0;250,128,114;255,255,0;173,255,47;0,0,255;255,0,255;0,255,255;210,180,140;255,255,224;0,0,0;64,224,208;0,128,0;25,25,112;224,255,255;165,42,42;255,192,203;128,0,128];

% compute correlations with gene expression:
mdl=[];
for mi = 1:19
    mdl{mi,1} = fitlm(MRIPar(:,1), Gene(:,mi));
    mdl{mi,2} = fitlm(MRIPar(:,2), Gene(:,mi));
    mdl{mi,3} = fitlm(MRIPar(:,3), Gene(:,mi));
    
    pval(mi,1)=mdl{mi,1}.Coefficients{2,4};
    pval(mi,2)=mdl{mi,2}.Coefficients{2,4};
    pval(mi,3)=mdl{mi,3}.Coefficients{2,4};
end
pALLtest=[pval];
[FDRval] = mafdr(pALLtest(:),'BHFDR', true);
FDRval=reshape(FDRval,size(pALLtest));

% find significant correlations:
signif{1}=find(FDRval(:,1)<0.05);
signif{2}=find(FDRval(:,2)<0.05);
signif{3}=find(FDRval(:,3)<0.05);


% plot significant correlations (SupFig. 11D-F):
set(0, 'DefaultFigureVisible', 'on');
CI=[];
Xd=[];
xlable_vec={'1st PC MDM','1st PC standard qMRI','MTV [fraction]'};
par_vec={'MDM','qMRI','MTV'};
for ii=1:length(signif)
    tmp_signif=signif{ii};
    for jj=1:length(tmp_signif)
        mi=tmp_signif(jj);
        % calculate CI:
        figure;
        h=plot(mdl{mi,ii});
        CI{jj,ii}(:,1)=h(4).YData;
        CI{jj,ii}(:,2)=h(3).YData;
        Xd{jj,ii}=h(3).XData;
        close(gcf)
        m=figure;
        hold on
        % plot fitted line:
        linvec=[min(MRIPar(:,ii)),max(MRIPar(:,ii))];
        linePre=linvec.*mdl{mi,ii}.Coefficients{2,1}+mdl{mi,ii}.Coefficients{1,1};
        [Ypred,p22] = predict(mdl{mi,ii},MRIPar(:,ii),'Prediction','curve' );
        linePre2=Xd{jj,ii}.*mdl{mi,ii}.Coefficients{2,1}+mdl{mi,ii}.Coefficients{1,1};
        [hl,hp]=boundedline(Xd{jj,ii},linePre2,[linePre2'-CI{jj,ii}(:,2),CI{jj,ii}(:,1)-linePre2'],'alpha','transparency',0.2,'nan','fill','cmap',ModuleColor(mi,:),'orientation','vert') ;
        delete(hl);
        % plot data points:
        scatter(MRIPar(:,ii),Gene(:,mi),20,ModuleColor(mi,:),'filled');
        grid on
        Module_name = modules{mi};
        plot(linvec,linePre,'--','Color', ModuleColor(mi,:), 'LineWidth',2)
        xlabel(xlable_vec{ii});
        ylabel([Module_name(1:end) ' Module' ]);
        % plot significance:
        [stars]=pval2stars(FDRval(mi,ii),'num');
        if ii==1
            xlim([-4.5 4]);
        elseif ii==2
            xlim([-4 4]);
        end
        set(gca,'FontSize',15);
        eq={strcat('R^2=',num2str(mdl{mi,ii}.Rsquared.Adjusted,'%0.2g')),stars};
        ypos=max(get(gca, 'ylim'));
        xpos=min(get(gca, 'xlim'));
        text(xpos-0.05*xpos,ypos-0.17*ypos,eq,'FontSize',14)
        hold off
        set(gcf, 'Position',[1 1 1268/2 460]);
    end
end

%% SupFig. 9B-C: The two gene modules most correlated with MDM measurements
figure;
x=1;
jj_tmp=[1,3];
for mi=[1,4]
    subplot(1,2,x)
    ii=1;
    hold on
    jj=jj_tmp(x);
    % plot confidence interval:
    linvec=[min(MRIPar(:,ii)),max(MRIPar(:,ii))];
    linePre=linvec.*mdl{mi,ii}.Coefficients{2,1}+mdl{mi,ii}.Coefficients{1,1};
    [Ypred,p22] = predict(mdl{mi,ii},MRIPar(:,ii),'Prediction','curve' );
    linePre2=Xd{jj,ii}.*mdl{mi,ii}.Coefficients{2,1}+mdl{mi,ii}.Coefficients{1,1};
    [hl,hp]=boundedline(Xd{jj,ii},linePre2,[linePre2'-CI{jj,ii}(:,2),CI{jj,ii}(:,1)-linePre2'],'alpha','transparency',0.2,'nan','fill','cmap',ModuleColor(mi,:),'orientation','vert') ;
    delete(hl);
    % plot data point:
    scatter(MRIPar(:,ii),Gene(:,mi),20,ModuleColor(mi,:),'filled');
    grid on
    Module_name = modules{mi};
    plot(linvec,linePre,'--','Color', ModuleColor(mi,:), 'LineWidth',2)
    xlabel(xlable_vec{ii});
    ylabel([Module_name(1:end) ' Module' ]);
    % plot significance:
    [stars]=pval2stars(FDRval(mi,ii),'num');
    t_leverage = 2*mdl{mi,ii}.NumCoefficients/mdl{mi,ii}.NumObservations;
    rm=find(mdl{mi,ii}.Diagnostics.Leverage > t_leverage);
    x1=MRIPar(:,ii);
    x1(rm)=[];
    y=Gene(:,mi);
    y(rm)=[];
    md2{x}=fitlm(x1,y);
    if ii==1
        xlim([-4.5 4]);
    elseif ii==2
        xlim([-4 4]);
    end
    set(gca,'FontSize',15);
    eq={strcat('R^2=',num2str(mdl{mi,ii}.Rsquared.Adjusted,'%0.2g')),stars};
    ypos=max(get(gca, 'ylim'));
    xpos=min(get(gca, 'xlim'));
    text(xpos-0.05*xpos,ypos-0.17*ypos,eq,'FontSize',15)
    hold off
    x=x+1;
end
set(gcf, 'Position',[1 1 1268*1.5*2/3 500]);

%%  SupFig 12 : The effect of outliers on the correlation of MDM with the gene moduls

figure;
% sup. fig. 12A: Identifying outliers based on high leverage - 1st module:
subplot(2,2,1)
mdl{1,1}.plotDiagnostics
xlabel('cortical ROIs');
set(gca,'FontSize',15);
xlim([0 65]);
title('');

% sup. fig. 12A: Identifying outliers based on high leverage - 2nd module:
subplot(2,2,2)
mdl{4,1}.plotDiagnostics
xlabel('cortical ROIs');
set(gca,'FontSize',15);
xlim([0 65]);
title('');

% sup. fig. 12BC: The regression between the 1st PC of MDM and two gene modules after removing the seven outliers identified in A.
x=1;
jj_tmp=[1,3];
for mi=[1,4]
    subplot(2,2,x+2)
    ii=1;
    hold on
    jj=jj_tmp(x);
    x1=MRIPar(:,ii);
    % find high laverage points:
    t_leverage = 2*mdl{mi,ii}.NumCoefficients/mdl{mi,ii}.NumObservations;
    % romove these points:
    rm=find(mdl{mi,ii}.Diagnostics.Leverage > t_leverage);
    x1(rm)=[];
    y=Gene(:,mi);
    y(rm)=[];
    % find CI:
    fig_tmp=figure;
    h=plot(md2{x});
    CI2(:,1)=h(4).YData;
    CI2(:,2)=h(3).YData;
    Xd2=h(3).XData;
    close(fig_tmp)
    % plot fitted line:
    linvec=[min(x1),max(x1)];
    linePre=linvec.*md2{x}.Coefficients{2,1}+md2{x}.Coefficients{1,1};
    [Ypred,p22] = predict(md2{x},x1,'Prediction','curve' );
    linePre2=Xd2.*md2{x}.Coefficients{2,1}+md2{x}.Coefficients{1,1};
    [hl,hp]=boundedline(Xd2,linePre2,[linePre2'-CI2(:,2),CI2(:,1)-linePre2'],'alpha','transparency',0.2,'nan','fill','cmap',ModuleColor(mi,:),'orientation','vert') ;
    delete(hl);
    % plot data points:
    scatter(x1,y,20,ModuleColor(mi,:),'filled');
    grid on
    Module_name = modules{mi};
    plot(linvec,linePre,'--','Color', ModuleColor(mi,:), 'LineWidth',2)
    xlabel(xlable_vec{ii});
    ylabel([Module_name(1:end) ' Module' ]);
    [stars]=pval2stars(md2{x}.Coefficients.pValue(2),'num');
    if ii==1
        xlim([-4.5 4]);
    elseif ii==2
        xlim([-4 4]);
    end
    set(gca,'FontSize',15);
    % plot significance:
    eq={strcat('R^2=',num2str(md2{x}.Rsquared.Adjusted,'%0.2g')),stars};
    ypos=max(get(gca, 'ylim'));
    xpos=min(get(gca, 'xlim'));
    text(xpos-0.05*xpos,ypos-0.17*ypos,eq,'FontSize',15)
    hold off
    x=x+1;
end
set(gcf, 'Position',[1 1 1268*1.5*2/3 500*2]);

checkBenferroni(1)=md2{1}.Coefficients.pValue(2)<0.05/(length(mdl(:))+3); 
checkBenferroni(2)=md2{2}.Coefficients.pValue(2)<0.05/(length(mdl(:))+3); 

end
