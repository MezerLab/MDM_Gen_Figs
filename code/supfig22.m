
function supfig22(ROI_list,v,fit_info,young_ind,l,comrl,r,color_new)

% this function generates sup. fig. 22: Distinct regional dependencies on MTV following R2* correction. 
vec=v;
ROI_list_new=ROI_list;
% plot figure
figure;
for j=1:4 % loop over qMRI parameters
    legend_info=[];
    % take young subjects data:
    slope_sub=nan(length(young_ind).*2,length(vec)); % subjects X ROIs
    for ii=1:length(vec) % loop over left ROIs 
        part=vec(ii);
        qMRIfit=fit_info.data{j};
        slope_sub(1:length(young_ind),ii)=squeeze(qMRIfit(part,1,young_ind));
        if ismember(part,l) && comrl % join right ROIs
            ind=find(l==part);
            slope_subr=squeeze(qMRIfit(r(ind),1,young_ind));
            slope_sub(length(young_ind)+1:2*length(young_ind),ii)=slope_subr;
        end
    end
    % use the median over subjects to arrange ROIs:
    med=nanmedian(slope_sub,1);
    subplot(3,2,(j))
    hold on
    [B,I]=sort(med);
    slope_sub=slope_sub(:,I);
    mat_color=color_new((I));
    newORd=vec;
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

subplot(3,2,6)
legend_info=[];
hold on
for i=1:length(v)
    plot([0,0],[0,0],'LineWidth',2,'color',color_new{i});
    legend_info{i}=(ROI_list{v((i))});
end
axis off
legend(legend_info)
set(get(gca,'legend'),'FontSize',15)
legend('boxoff')
set(gcf, 'Position',[1 1 1030 1500]);
legend(legend_info,'Location','bestoutside');

end
