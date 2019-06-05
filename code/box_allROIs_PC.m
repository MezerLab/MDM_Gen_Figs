
function box_allROIs_PC(hip,color_new,FDRval,ROI_list)

% this function generates SupFig. 18:   aging related changes revealed by the 1st principle component of MDM. 

coef=[1,3];
hipPC=[];
j=1;

% comput PCA
for par=1:length(coef) % MDM,qMRI
    for age=1:2
        MDMmat=nanmedian(squeeze((hip(:,age,coef(par),:,:))),1);
        MDMmat=squeeze(MDMmat);
        [coeffMDMyoung,scoreMDMyoung] = pca(zscore(MDMmat));
        for ii=1:size(hip,1)
            hipPC(ii,age,coef(par),:,:)=squeeze((hip(ii,age,coef(par),:,:)))*coeffMDMyoung(:,1);
        end
    end

end
hipPC(:,:,2,:)=hip(:,:,2,:,2);

for par=1:3 % MDM,qMRI
for roi_ind=1:size(hip,4)
    [hALLPC(par,roi_ind),pALLPC(par,roi_ind)] = ttest2(hipPC(:,1,par,roi_ind),hipPC(:,2,par,roi_ind));
end
end

% correct for multiple comparisions using FDR
pALLtestPC=pALLPC;
[FDRvalPC] = mafdr(pALLtestPC(:),'BHFDR', true);
FDRvalPC=reshape(FDRvalPC,size(pALLPC));
HvalPC=FDRvalPC;
HvalPC(FDRvalPC>(5./100))=nan;

% plot figure
tot=[7,14];
tot2=[1,8];
ind_subplot(:,1)=[1:3,7:9,13:15,19:21,25:27,31:33,37:39]';
ind_subplot(:,2)=[4:6,10:12,16:18,22:24,28:30,34:36,40:42];
h=figure;
for k=1:length(tot)
    sub_n=1;
    for ii=tot2(k):tot(k)
        x=1;
        subplot(8,6,ind_subplot(sub_n,k))
        tmp=reshape(hipPC(:,:,1,ii),[size(hipPC(:,:,1,ii),1),1,size(hipPC(:,:,1,ii),2)]);
        m=iosr.statistics.boxPlot(tmp);
        m.boxColor={color_new{(ii)},  [ 0.7529  0.7529  0.7529]};
        m.boxAlpha=0.5;
        m.medianColor='k';
        m.lineColor={color_new{(ii)}, [ 0.7529  0.7529  0.7529]};
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
        if x==1 && FDRvalPC(1,ii)<=(5./100)
            pvalm=FDRvalPC(1,ii);
            if pvalm<=1E-3
                stars='***';
            elseif pvalm<=1E-2
                stars='**';
            elseif pvalm<=0.05
                stars='*';
            end
            text(0.45,0.95,[stars],'Units' ,'normalized')
        end
        if ind_subplot(sub_n,k)==1 || ind_subplot(sub_n,k)==4
            title(['1st PC of MDM']);
        end
        x=x+1;
        sub_n=sub_n+1;
        subplot(8,6,ind_subplot(sub_n,k))
        hold on
        tmp=reshape(hipPC(:,:,2,ii),[size(hipPC(:,:,2,ii),1),1,size(hipPC(:,:,2,ii),2)]);
        m=iosr.statistics.boxPlot(tmp);
        m.boxColor={color_new{(ii)},  [ 0.7529  0.7529  0.7529]};
        m.boxAlpha=0.5;
        m.medianColor='k';
        m.lineColor={color_new{(ii)}, [ 0.7529  0.7529  0.7529]};
        m.medianColor='k';
        m.notchLineColor='k';
        m.lineWidth=1.5;
        m.notch=true;
        m.limit='none';
        m.scatterMarker='X';
        m.scatterSize=10;
        m.showScatter=false;
        xticklabels=[];
        m.groupWidth=1.5;
        set(gca, 'xtick',[0.625 1.375], 'xticklabel', xticklabels);
        set(gcf,'color','w');
        if x==2 && FDRval(2,ii,2)<=(5./100)
            pvalm=FDRval(2,ii,2);
            if pvalm<=1E-3
                stars='***';
            elseif pvalm<=1E-2
                stars='**';
            elseif pvalm<=0.05
                stars='*';
            end
            text(0.45,0.95,[stars],'Units' ,'normalized')
        elseif x==2 && FDRval(2,ii,2)>(5./100)
        end
        if ind_subplot(sub_n,k)==2 || ind_subplot(sub_n,k)==5
            title(['MTV']);
        end
        x=x+1;
        
        sub_n=sub_n+1;
        subplot(8,6,ind_subplot(sub_n,k))
        tmp=reshape(hipPC(:,:,3,ii),[size(hipPC(:,:,3,ii),1),1,size(hipPC(:,:,3,ii),2)]);
        m=iosr.statistics.boxPlot(tmp);
        m.boxColor={color_new{(ii)},  [ 0.7529  0.7529  0.7529]};
        m.boxAlpha=0.5;
        m.medianColor='k';
        m.lineColor={color_new{(ii)}, [ 0.7529  0.7529  0.7529]};
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
        if x==3 && FDRvalPC(3,ii)<=(5./100)
            pvalm=FDRvalPC(3,ii);
            if pvalm<=1E-3
                stars='***';
            elseif pvalm<=1E-2
                stars='**';
            elseif pvalm<=0.05
                stars='*';
            end
            text(0.45,0.95,[stars],'Units' ,'normalized')
        elseif x==3 && FDRvalPC(3,ii)>(5./100)
        end
        if ind_subplot(sub_n,k)==3 || ind_subplot(sub_n,k)==6
            title('1st PC of qMRI');
        end
        sub_n=sub_n+1;
        
        
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