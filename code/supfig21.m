
function supfig21(R2mat,color_new,ROI_list,v,FDRval)

% this function generate sup. fig. 21: Iron-related changes between the two age groups.

% p-values corrected for multiple comparisions:
FDR_R2s=FDRval(:,:,5);

% plot figure
h=figure;
for ii=1:size(R2mat,4)
    subplot(4,4,ii)
    tmp=reshape(R2mat(:,:,3,(ii)),[size(R2mat(:,:,3,(ii)),1),1,size(R2mat(:,:,3,(ii)),2)]);
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
    fdr_tmp=FDR_R2s(3,ii);
    if  fdr_tmp<=(5./100)
        [stars]=pval2stars(fdr_tmp,'stars');
        text(0.35,0.95,[stars],'Units' ,'normalized')
    end
    ylabel('R2* s^-^1]');
    title(ROI_list(v(ii)),'FontWeight','Normal');
    set(gca,'FontSize',15);

end
set(gcf, 'Position',[1 1 1250 1050]);

end

