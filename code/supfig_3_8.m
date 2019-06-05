function supfig_3_8(main_dir,ROI_list,fit_info,v,l,r,color_new,comrl,young_ind)

% this function generates sup. figures 3-8

%% SupFig. 3:  Linear dependency of R1 and MTsat on MTV in the human brain.


% load single subject data
load(fullfile(main_dir,'/fig2_5_6/human_data_2.mat'));
dist=[distMT ; distR1];
% dist is a struct containing the data used for the fit for each of the 35 ROIs and for each qMRI
% parameter (R1 and MTsat).

fit=fitMT;
fit(:,:,2)=fitR1;
% fit is a 35X6X2 matrices. There are 35 ROIs and 2 qMRI parameters (R1 and MTsat). The 2nd dimention represent: [slope,
% intersection, mean MTV, mean qMRI parameter, R^2 and cvRMSE] for each ROI.


roi=v;
range={[0 0.045]};
for ii=1:size(fit,3) % loop over qMRI parameters
    h=figure('Position', get(0, 'Screensize'));
    for i=1:length(roi) % loop over ROIs
        subplot(4,4,i)
        hold on
        v1=dist(ii,roi(i)).xbin;
        v2=dist(ii,roi(i)).ybin;
        % plot data points
        plot(v1,v2,'.','color',color_new{(i)},'MarkerSize',15)
        err=dist(ii,roi(i)).err;
        % plot error for each bin:
        [hl,hp]=boundedline(v1,v2,err,'.','alpha','transparency',0.2,'cmap',color_new{(i)}) ;
        % plot fitted line:
        plot(v1,v1*fit(roi(i),1,ii) +fit(roi(i),2,ii),'color',color_new{(i)},'LineWidth',2);
        hold off
        xlabel('MTV [fraction]')
        if ii==2
            ylabel('R1 [s^-^1]')
        elseif ii==1
            ylabel('MTsat [p.u.]')
        elseif ii==3
            ylabel('MD [μm^2/ms]')
        elseif ii==4
            ylabel('R2 [s^-^1]')
        end
    end
end
%% SupFig. 5 - Variance explained by the qMRI dependencies on MTV and its stability between hemispheres.

% sup fig 5A: Distributions of adjusted R^2
h=figure;
for ii=1:length(fit_info.str) % loop over qMRI parameters
    subplot(2,2,ii)
    data_m=fit_info.data{ii};
    data_m=squeeze(data_m(:,5,:));
    histogram(data_m(:)); % find R^2 for this parameter
    title(['R^2 for d' fit_info.str{ii} ' /dMTV' ]);
    eq=['median R^2=' num2str(nanmedian(data_m(:)),'%.2g') '±' num2str(mad(data_m(:)),'%.2g')];
    disp(nanmedian(data_m(:)));
    annotation('textbox',[0.166 0.72 0.9 0.15],'EdgeColor','none','String',eq,'FontSize',15)
    set(gca,'FontSize',15);
    
end


%  sup fig 5B: Agreement in the MTV derivatives of qMRI parameters between hemispheres
h=figure;
for ii=1:length(fit_info.str)
    subplot(2,2,ii)
    data_m=fit_info.data{ii};
    data_m=squeeze(data_m(:,1,:));
    Rh=data_m(r,:);
    Lh=data_m(l,:);
    mdl = fitlm(Rh(:),Lh(:));
    R=mdl.Rsquared.Adjusted;
    P=mdl.Coefficients{2,4};
    hold on
    plot(Rh(:),Lh(:),'*')
    hline = refline(1,0);
    hline.Color = 'k';
    hline.LineStyle='--';
    xlabel(['rh d' fit_info.str{ii} '/dMTV' ]);
    ylabel(['lh d' fit_info.str{ii} ' /dMTV' ]);
    
    sl=mdl.Coefficients{2,1};
    int=mdl.Coefficients{1,1};
    [stars]=pval2stars(mdl.Coefficients{2,4},'num');
    eq=['R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f'),' ',[stars]];
    title(eq)
    set(gca,'FontSize',15);
end
%% SupFig. 6- Distinct regional dependencies on MTV for R2 and MD.

vec=v;
ROI_list_new=ROI_list;
s=figure;
sub_ind=[ 0 0 1 2];
for j=3:4
    legend_info=[];
    % take young subjects' MDM measurements:
    slope_sub=nan(length(young_ind).*2,length(vec));
    for ii=1:length(vec) % for left ROIs
        part=vec(ii);
        qMRIfit=fit_info.data{j};
        slope_sub(1:length(young_ind),ii)=squeeze(qMRIfit(part,1,young_ind));
        if ismember(part,l) && comrl % and right ROIs
            ind=find(l==part);
            slope_subr=squeeze(qMRIfit(r(ind),1,young_ind));
            slope_sub(length(young_ind)+1:2*length(young_ind),ii)=slope_subr;
        end
    end
    % avarage over subjects:
    med=nanmedian(slope_sub,1);
    subplot(1,2,sub_ind(j))
    hold on
    [B,I]=sort(med); % sort ROIs by median:
    slope_sub=slope_sub(:,I);
    mat_color=color_new((I));
    ROI_list_new_sort=ROI_list_new((I));
    newORd=vec;
    % plot figure:
    for i=1:length(vec)
        plot([0,0],[0,0],'LineWidth',2,'color',color_new{i});
        legend_info{i}=(ROI_list_new{(newORd((i)))});
    end
    m=[];
    m=iosr.statistics.boxPlot(reshape(slope_sub,[length(young_ind).*2,1,length(vec)]));
    m.boxColor=mat_color;
    m.groupWidth=1;
    m.lineColor=mat_color;
    m.medianColor='k';
    m.notch=false;
    m.notchLineColor='k';
    m.lineWidth=1;
    m.limit='none';
    grid on
    if j==1
        ylabel(['d' fit_info.str{j} '/dMTV [p.u.]']);
    elseif j==3
        ylabel(['d' fit_info.str{j} '/dMTV [μm^2/ms]']);
    else
        ylabel(['d' fit_info.str{j} '/dMTV [S^-^1]']);
    end
    hold off
    set(gca, 'xtick', [], 'xticklabel', []);
    set(gcf,'color','w');
    set(gca,'FontSize',15)
    
end
set(gcf, 'Position',[1 1 1030 500]);
legend(legend_info,'Location','bestoutside');



%% SupFig. 8 -  Standard qMRI parameters’ variation across the brain. 

vec=v;
ROI_list_new=ROI_list;
s=figure;
sub_ind=[ 0 0 1 2];
for j=1:length(fit_info.str) % loop over qMRI parameters
    legend_info=[];
    slope_sub=nan(length(young_ind).*2,length(vec)); % take young subjects data
    for ii=1:length(vec) % loop over left ROIs
        part=vec(ii);
        qMRIfit=fit_info.data{j};
        % find mean value of qMRI parameters in each ROI:
        slope_sub(1:length(young_ind),ii)=squeeze(qMRIfit(part,4,young_ind));
        if ismember(part,l) && comrl % join right ROIs
            ind=find(l==part);
            slope_subr=squeeze(qMRIfit(r(ind),4,young_ind));
            slope_sub(length(young_ind)+1:2*length(young_ind),ii)=slope_subr;
        end
    end
    % average over subjects:
    med=nanmedian(slope_sub,1);
    % plot figure:
    subplot(2,2,(j))
    hold on
    % sort ROIs based on median values:
    [B,I]=sort(med);
    slope_sub=slope_sub(:,I);
    mat_color=color_new((I));
    ROI_list_new_sort=ROI_list_new((I));
    newORd=[1:5,8:12,7,6,13:19];
    newORd=vec;
    % prepare legend:
    for i=1:length(vec)
        plot([0,0],[0,0],'LineWidth',2,'color',color_new{i});
        legend_info{i}=(ROI_list_new{(newORd((i)))});
    end
    % plot boxplot:
    m=[];
    m=iosr.statistics.boxPlot(reshape(slope_sub,[length(young_ind).*2,1,length(vec)]));
    m.boxColor=mat_color;
    m.groupWidth=1;
    m.lineColor=mat_color;
    m.medianColor='k';
    m.notch=false;
    m.notchLineColor='k';
    m.lineWidth=1;
    m.limit='none';
    grid on
    if j==1
        ylabel([fit_info.str{j} ' [p.u.]']);
    elseif j==3
        ylabel([fit_info.str{j} ' [μm2/ms]']);
    else
        ylabel([fit_info.str{j} ' [S^-^1]']);
    end
    hold off
    set(gca, 'xtick', [], 'xticklabel', []);
    set(gcf,'color','w');
    set(gca,'FontSize',15)
end
set(gcf, 'Position',[1 1 1030 1000]);
legend(legend_info,'Location','bestoutside');

end
