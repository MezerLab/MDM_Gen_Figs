
function inverseMTV(ROI_list,fit_info,v,l,r,color_new,comrl,young_ind)

% this function generates sup. fig. 30: Examination of the linear relationship of qMRI parameters and 1/PD. 

% ROI to use:
vec=v;
ROI_list_new=ROI_list;

figure;
ROI_list_new{1}='Thalamus';
ROI_list_new{4}='Pallidum';
ROI_list_new{13}='WM-Frontal';


%% Sup. Fig. 30C:  Similarly to MTV, the 1/WC derivatives of qMRI parameters are region-specific

for j=1:2 % loop over qMRI parameters
    subplot(1,2,j)
    legend_info=[];
    slope_sub=nan(length(young_ind).*2,length(vec));
    for ii=1:length(vec) % loop over ROI
        part=vec(ii);
        qMRIfit=fit_info.data{j};
        % take measurements of young subjects in the left ROI
        slope_sub(1:length(young_ind),ii)=squeeze(qMRIfit(part,1,young_ind));
        if ismember(part,l) && comrl % add right ROIs
            ind=find(l==part);
            slope_subr=squeeze(qMRIfit(r(ind),1,young_ind));
            slope_sub(length(young_ind)+1:2*length(young_ind),ii)=slope_subr;
        end
    end
    % use the median over each ROI to determine the order:
    med=nanmedian(slope_sub,1);
    hold on
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
    % plot boxplot
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
        ylabel(['d' fit_info.str{j} '/d(1/WC) [p.u.]']);
        ylim([-0.1 11]);
        set(gca, 'ytick', [0 2.5 5 7.5 10]);
    elseif j==3
        ylabel(['d' fit_info.str{j} '/d(1/WC) [μm^2/ms]']);
    else
        ylabel(['d' fit_info.str{j} '/d(1/WC) [S^-^1]']);
    end
    hold off
    set(gca, 'xtick', [], 'xticklabel', []);
    set(gcf,'color','w');
    set(gca,'FontSize',15)
end
set(gcf, 'Position',[1 1 1040 970/2]);

%% Sup. Fig. 30D: In a biological range, 1/WC varies almost linearly with MTV. 

% biological range of MTV values:
MTV=linspace(0.1,0.4,100);

figure
subplot(1,2,1)
plot(MTV,1./(1-MTV),'-.','LineWidth',3)
MTV=linspace(0.00001,1,100);
xlim([0.095,0.41]);
grid on
xlabel('MTV [fraction]');
ylabel('1/WC [fraction]');
set(gca,'FontSize',15);
subplot(1,2,2)
plot(MTV,1./(1-MTV),'-.','LineWidth',3)
grid on
set(gcf, 'Position',[1 1 1040 970./2]);
xlim([-0.1,1.05]);
ylim([-3,103]);


end