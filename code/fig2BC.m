
function fig2BC(main_dir,ROI_list,fit_info,v,l,r,color_new,comrl,young_ind)

% This funcition generates figure 2BC

% Fig. 2B : Calculating the MDM signatures. 

vec=v;
ROI_list_new=ROI_list;

% load single subject data
load(fullfile(main_dir,'/fig2_5_6/human_data_2.mat'));
dist=[distMT ; distR1];
% dist is a struct containing the data used for the fit for each of the 35 ROIs and for each qMRI
% parameter (R1 and MTsat).

fit=fitMT;
fit(:,:,2)=fitR1;
% fit is a 35X6X2 matrices. There are 35 ROIs and 2 qMRI parameters (R1 and MTsat). The 2nd dimention represent: [slope,
%   intersection, mean MTV, mean qMRI parameter, R^2 and cvRMSE] for each ROI.


% plot figure:
figure;
ROI_list_new{1}='Thalamus';
ROI_list_new{4}='Pallidum';
ROI_list_new{13}='WM-Frontal';

for ii=[1,2]
    subplot(2,2,ii)
    roi=[1,11,4];
    legendInfo=[];
    hold on
    % prepare legend:
    for i=1:length(roi)
        plot([0,0],[0,0],'LineWidth',2,'color',color_new{roi(i)});
        legendInfo{i}=strcat(ROI_list_new{vec(roi(i))});
    end
    % plot data points:
    for i=1:length(roi)
        v1=dist(ii,vec(roi(i))).xbin;
        v2=dist(ii,vec(roi(i))).ybin;
        plot(v1,v2,'.','color',color_new{roi(i)},'MarkerSize',15)
    end
    % plot the  median absolute deviation for each bin:
    for i=1:length(roi)
        v1=dist(ii,vec(roi(i))).xbin;
        v2=dist(ii,vec(roi(i))).ybin;
        err=dist(ii,vec(roi(i))).err;
        [hl,hp]=boundedline(v1,v2,err,'.','alpha','transparency',0.2,'cmap',color_new{roi(i)}) ;
        
    end
    % plot fitted line:
    for i=1:length(roi)
        v1=dist(ii,vec(roi(i))).xbin;
        plot(v1,v1*fit(vec(roi(i)),1,ii) +fit(vec(roi(i)),2,ii),'color',color_new{(roi(i))},'LineWidth',2);
    end
    grid on
    hold off
    xlabel('MTV [fraction]')
    set(gca,'FontSize',15)
    if ii==2
        legend(legendInfo,'Location','northwest')
        xlim([0.2135 0.37])
        ylim([0.765 1.27])
        ylabel('R1 [S^-1]')
    else
        xlim([0.2135 0.37])
        ylim([0.034 0.059]);
        ylabel('MTsat [p.u.]')
        tick=get(gca, 'ytick');
        set(gca, 'ytick', tick, 'yticklabel', tick*100);
    end
end

%% Fig. 2C:  The reliability of the MDM method across subjects. 

% take only young subjects' data and join estimates for left and right
% ROIs:
for j=1:2
    subplot(2,2,j+2)
    legend_info=[];
    slope_sub=nan(length(young_ind).*2,length(vec));
    for ii=1:length(vec) % loop over ROIs
        part=vec(ii);
        % take MDM measuremet for left ROIs:
        qMRIfit=fit_info.data{j}; 
        slope_sub(1:length(young_ind),ii)=squeeze(qMRIfit(part,1,young_ind));
        % join right ROIs as well:
        if ismember(part,l) && comrl
            ind=find(l==part);
            slope_subr=squeeze(qMRIfit(r(ind),1,young_ind));
            slope_sub(length(young_ind)+1:2*length(young_ind),ii)=slope_subr;
        end
    end
    
    % take the median over all subjects:
    med=nanmedian(slope_sub,1);
    hold on
    
    % ararnge ROIs by their median values:
    [B,I]=sort(med);
    slope_sub=slope_sub(:,I);
    mat_color=color_new((I));
    newORd=vec;
    
    % boxplot:
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
    pos = get(gca, 'Position');
    pos(2)=0.15;
    set(gca, 'Position', pos)
    if j==1
        ylabel(['d' fit_info.str{j} '/dMTV [p.u.]']);
        ylim([-0.1 21]);
        set(gca, 'ytick', [0 5 10 15 20]);
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
set(gcf, 'Position',[1 1 1040 970]);
%% Fig 2C:  legend of ROIs

ROI_list=strrep(ROI_list,'ctx','CTX');
ROI_list=strrep(ROI_list,'wm','WM');
figure
legend_info=[];
hold on
for i=1:length(v)
    plot([0,0],[0,0],'LineWidth',2,'color',color_new{i});
    legend_info{i}=(ROI_list{v((i))});
end
legend(legend_info)
set(get(gca,'legend'),'FontSize',12)
legend('boxoff')
end