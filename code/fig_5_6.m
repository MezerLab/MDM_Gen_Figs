

function fig_5_6(hip,color_new,Test_young,Test_old,FDRval)

% this function generates Fig. 5-6: Region-specific aging-related molecular changes revealed by the R1/MTsat dependency on MTV. 

vecALL=[7 5 13  ; 12, 1, 8]; % ROIs to show
for j=1:2 % loop over qMRI parameters
    h=figure;
    set(gcf, 'Position',[1 1 1200 1050]);
    x=1;
    vec=vecALL(j,:);
    for ii=1:length(vec) % loop over ROIs
        for par=[3,1,2]
            % plot boxplot:
            h(par)=subplot(3,4,x);
            tmp=reshape(hip(:,:,par,vec(ii),j),[size(hip(:,:,par,vec(ii),j),1),1,size(hip(:,:,par,vec(ii),j),2)]);
            m=iosr.statistics.boxPlot(tmp);
            m.boxColor={color_new{(vec(ii))},  [ 0.7529  0.7529  0.7529]};
            m.boxAlpha=0.5;
            m.medianColor='k';
            m.lineColor={color_new{(vec(ii))}, [ 0.7529  0.7529  0.7529]};
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
            NumTicks = 3;
            L = get(gca,'YLim');
            set(gca,'YTick',linspace(L(1),L(2),NumTicks))
            if par==2
                fdr_tmp=FDRval(par,vec(ii),2);
            else
                fdr_tmp=FDRval(par,vec(ii),j);
            end
            if  fdr_tmp<=(5./100)
                [stars]=pval2stars(fdr_tmp,'stars');
                text(0.35,0.95,[stars],'Units' ,'normalized')
            end
            x=x+1;
            if par==3
                if j==1
                    title(['MTsat [p.u.]']);
                elseif j==2
                    title(['R1 [s^-^1]']);
                end
            elseif par==1
                if j==1
                    title(['dMTsat/dMTV [p.u.]']);
                elseif j==2
                    title(['dR1/dMTV [s^-^1]']);
                end
            elseif par==2
                title(['MTV [fraction]']);
            end
        end
        
        set(findall(h, 'Type', 'Text'),'FontWeight', 'Normal')
        
        % spider plot:
        color_spider=color_new;
        h(4)=subplot(3,4,x);
        curr=[];
        opt_axes=[];
        opt_lines=[];
        opt_area=[];
        curr(:,:,1)=squeeze(Test_young(:,vec(ii),:));
        curr(:,:,2)=squeeze(Test_old(:,vec(ii),:));
        margins= [curr(:,:,1)  curr(:,:,2)];
        curr=permute(curr,[1 3 2]);
        Multi(:,1)= prctile(margins',5)';
        Multi(:,2)= prctile(margins',95)';
        opt_area.err        = 'std';
        opt_area.FaceAlpha  = 0.5;
        opt_area.Color      = [color_spider{vec(ii)}; [175,175,175]./255];
        opt_lines.Color     = [color_spider{vec(ii)}; [175,175,175]./255];
        opt_lines.LineWidth = 2;
        opt_lines.LineStyle = '-';
        opt_lines.Marker    = 'none';
        opt_lines.Labels    = false;
        opt_lines.Legend    = {'younger','older'};
        opt_axes.Signif=squeeze(FDRval(1,vec(ii),[1 2 3 4]));
        opt_axes.Background = 'w';
        opt_axes.Multi   = Multi;
        opt_axes.Labels     = {'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'};
        polygonplotMDM(curr,opt_axes,opt_lines,opt_area);
        set(gca,'FontSize',15)
        x=x+1;
        pos1=get(h(3),'Position');
        pos1(3)=0.7*pos1(3);
        pos2=get(h(1),'Position');
        pos2(1)=pos2(1)-0.12*pos2(1);
        pos2(3)=0.7*pos2(3);
        pos3=get(h(2),'Position');
        pos3(1)=pos3(1)-0.18*pos3(1);
        pos3(3)=0.7*pos3(3);
        set(h(3), 'Position', pos1);
        set(h(2), 'Position', pos3);
        set(h(1), 'Position', pos2);
        pos4=get(h(4),'Position');
        pos4(3)=6*pos4(3);
        pos4(1)=pos4(1)-0.68*pos4(1);
        set(h(4), 'Position', pos4);
    end
end
end