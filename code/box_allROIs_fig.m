
function box_allROIs_fig(str_vec,hip,fit_info,ROI_list,color_new,FDRval)

% this function generates SupFig. 14-17: separating molecular and water related contributions. 

for j=1:length(str_vec) % loop over qMRI parameters
    tot=[7,14];
    tot2=[1,8];
    ind_subplot(:,1)=[1:3,7:9,13:15,19:21,25:27,31:33,37:39]';
    ind_subplot(:,2)=[4:6,10:12,16:18,22:24,28:30,34:36,40:42];
    h=figure;
    for k=1:length(tot)
        sub_n=1;
        for ii=tot2(k):tot(k)
            x=1;
            for par=[1,2,3]
                subplot(8,6,ind_subplot(sub_n,k))
                tmp=reshape(hip(:,:,par,(ii),j),[size(hip(:,:,par,(ii),j),1),1,size(hip(:,:,par,(ii),j),2)]);
                m=iosr.statistics.boxPlot(tmp);
                m.boxColor={color_new{((ii))},  [ 0.7529  0.7529  0.7529]};
                m.boxAlpha=0.5;
                m.medianColor='k';
                m.lineColor={color_new{((ii))}, [ 0.7529  0.7529  0.7529]};
                m.medianColor='k';
                m.notchLineColor='k';
                m.lineWidth=1.2;
                m.notch=true;
                m.limit='none';
                m.scatterMarker='X';
                m.scatterSize=10;
                m.showScatter=false;
                xticklabels=[];
                m.groupWidth=1.5;
                set(gca, 'xtick',[0.625 1.375], 'xticklabel', xticklabels);
                set(gcf,'color','w');
                if par==2
                    fdr_tmp=FDRval(par,(ii),2);
                else
                    fdr_tmp=FDRval(par,(ii),j);
                end
                if  fdr_tmp<=(5./100)
                    [stars]=pval2stars(fdr_tmp,'stars');
                    text(0.35,0.95,[stars],'Units' ,'normalized')
                end
                if ind_subplot(sub_n,k)==2 || ind_subplot(sub_n,k)==5
                    title('MTV');
                end
                if ind_subplot(sub_n,k)==3 || ind_subplot(sub_n,k)==6
                    title([fit_info.str{str_vec(j)}]);
                end
                if ind_subplot(sub_n,k)==1 || ind_subplot(sub_n,k)==4
                    title(['d' fit_info.str{str_vec(j)} '/dMTV']);
                end
                sub_n=sub_n+1;
                x=x+1;
            end
        end
    end
    
    subplot(8,6,46)
    legend_info=[];
    hold on
    vec1=vec(1:7);
    for i=1:length(vec1)
        plot([0,0],[0,0],'LineWidth',2,'color',color_new{i});
        legend_info{i}=(ROI_list{vec1((i))});
    end
    axis off
    legend(legend_info)
    set(get(gca,'legend'),'FontSize',12)
    legend('boxoff')
    
    subplot(8,6,48)
    legend_info=[];
    hold on
    vec2=vec(8:14);
    for i=1:length(vec2)
        plot([0,0],[0,0],'LineWidth',2,'color',color_new{7+i});
        legend_info{i}=(ROI_list{vec2((i))});
    end
    axis off
    legend(legend_info)
    set(get(gca,'legend'),'FontSize',12)
    legend('boxoff')
    set(gcf, 'Position',[1 1 1050 1050]);
    
end
end