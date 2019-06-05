function supfig1D(main_dir)

% This function generates supfig1D: The dependency of R1, R2 and MTsat on MTV for non-lipid compound  

%% bring data

Par={'R1 [s^-^1]','R2 [s^-^1]','MTsat [p.u.]'};
lip={'Sucrose','Glucose','BSA','BSA+Iron'};


load(fullfile(main_dir,'/fig1/phantom_data5.mat'));
% This .mat file contains:
% - data is a matrix of 4X3 : the rows are different compunds (see variable
%   lip) and the columns are different qMRI parameters (see variable Par).
% - MTV is a cell with 4 vectors corresponding to 4 different compounds (see variable
%   lip). In each vectore you can find the MTV values of the compound. 

%% Compute the linear fit between MTV and qMRI parameters for each compound:

fitobjectT=[];
for ii=1:length(lip)
    for jj=1:length(Par)
        [v1,I]=sort(MTV{ii});
        v2=data{ii,jj};
        v2=v2(I);
        v2=v2(~isnan(v1));
        v1=v1(~isnan(v1));
        fitobjectT{ii,jj} = fit(v1',v2','poly1');
        mdl{ii,jj} = fitlm(v1,v2);
        Rs(ii,jj)=mdl{ii,jj}.Rsquared.Adjusted;
    end
end

%% create figure 

colors=[252,141,89;...
    223,101,176;...
    127,191,123;...
    94,60,153]./255;

h=figure;
for ii=1:length(Par)
    subplot(1,3,ii)
    hold on
    for jj=1:length(lip)
        plot([0,0],[0,0],'LineWidth',2,'color',colors(jj,:));
        legendInfo{jj}=strcat(lip{jj});
    end
    for jj=1:length(lip)
        v1=MTV{jj};
        v2=data{jj,ii};
        plot(v1,v2,'.','color',colors(jj,:),'MarkerSize',15)
    end
    for jj=1:length(lip)
        [v1,I]=sort(MTV{jj});
        v2=data{jj,ii};
        v2=v2(I);
        v2=v2(~isnan(v1));
        v1=v1(~isnan(v1));
        fitobject = fitobjectT{jj,ii};
        p22 = predint(fitobject,v1',0.95,'functional','on');
        lineV=v1*fitobject.p1+fitobject.p2;
        [hl,hp]=boundedline(v1',lineV',[lineV'-p22(:,1),p22(:,2)-lineV'],'alpha','transparency',0.2,'nan','gap','cmap',colors(jj,:),'orientation','vert') ;
    end
    for jj=1:length(lip)
        fitobject = fitobjectT{jj,ii};
        [v1,I]=sort(MTV{jj});
        v1=v1(~isnan(v1));
        plot([0 0.4],[0 0.4]*fitobject.p1 +fitobject.p2,'--','color',colors((jj),:),'LineWidth',1);
        plot(v1,v1*fitobject.p1 +fitobject.p2,'color',colors((jj),:),'LineWidth',2);
    end
    hold off
    grid on
    xlim([-0.01 0.28])
    if ii==3
        ylim([-0.2 1.7]);
    elseif ii==1
        ylim([0.2 1.2]);
    elseif ii==2
        ylim([0 25]);
    end
    xlabel('MTV [fraction]')
    ylabel(Par{ii})
    set(gca,'FontSize',15);
end
legend(legendInfo,'Location','northwest')
set(gcf, 'Position',[1 1 1453 418]);

end 

