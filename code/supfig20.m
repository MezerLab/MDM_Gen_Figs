function supfig20(diff)

% this function generates sup. figure 20:  Distinct molecular aging trajectories- comparison of different MDM dimensions. 

% set color
c2=[149 143 83;149 170 146;  226 211 112; 222 145 103]./255;
c=[156 163 169; 161 162 120 ; 194 172 124; 152 136 146]./255;

%% supfig 20B: The spatial correlation between aging-related changes in different MTV derivatives. 

% compute correlations between aging of different MDM dimentions:
mdlL=[];
pvals=[];

% order of comparision between markers:
xf=[1,1,2,2,1,3];
yf=[3,4,3,4,2,4];
c3=[c(1,:) ; c2(2:4,:)];
for ii=1:6
    x=diff(1,:,xf(ii));
    y=diff(1,:,yf(ii));
    mdl = fitlm(x,y);
    % find CI:
    cifig=figure;
    h=plot(mdl);
    CI=[];
    CI(:,1)=h(4).YData;
    CI(:,2)=h(3).YData;
    Xd=h(3).XData;
    close(cifig)
    % generate prediction: 
    linvec=[min(x),max(x)];
    linePre=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1}; % fitted line
    linePre2=Xd.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1}; % prediction for each data points
    % save fit info:
    boundL{ii}=[linePre2'-CI(:,2),CI(:,1)-linePre2'];
    XdL{ii}=Xd;
    linePre2L{ii}=linePre2;
    linePreL{ii}=linePre;
    linvecL{ii}=linvec;
    mdlL{ii}=mdl;
    pvals(ii)=mdlL{ii}.Coefficients{2,4};
end

% plot figure:
sub=1;
lables={'dMTsat/dMTV Effect size','dR1/dMTV Effect size','dMD/dMTV Effect size','dR2/dMTV Effect size'};
figure
for ii=1:6
    subplot(3,2,ii)
    hold on
    x=diff(1,:,xf(ii));
    y=diff(1,:,yf(ii));
    [hl,hp]=boundedline(XdL{ii},linePre2L{ii},boundL{ii},'alpha','transparency',0.2,'nan','fill','cmap',[0.5 0.5 0.5],'orientation','vert') ;
    delete(hl);
    scatter(x,y,20,[0.5 0.5 0.5],'filled');
    plot(linvecL{ii},linePreL{ii},'--', 'LineWidth',2,'Color',[0.5 0.5 0.5])
    [stars]=pval2stars(pvals(ii),'num');
    eq={strcat('R^2=',num2str(mdlL{ii}.Rsquared.Adjusted,'%0.2g'),', ', stars)};
    xlabel(lables{xf(ii)});
    ylabel(lables{yf(ii)});
    title(eq)
    grid on
    sub=sub+1;
    hold off
    ax = gca;
    ax.XColor=c3(xf(ii),:);
    ax.YColor=c3(yf(ii),:);
    ax.FontWeight='bold';
end

set(gcf, 'Position',[1 1 1053 800]);

%% supfig 20A : Distinct spatial patterns of different molecular aging markers throughout the brain. 

% set order of ROIs
ord=[1:4,6,5,9,7,8,10,12,11,13,14];
[M,I] = max(abs(diff(1,:,:)),[],3);
lables={'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'};

% plot figure:
figure
hold on
for ii=1:14
    rectangle('Position',[(ii)-0.5 -2.5 1 8],'FaceColor',[c3(I(ii),:) 0.3],'EdgeColor','none')
end
rectangle('Position',[1-0.75 -0.2 14.76 0.4],'FaceColor',[1 1 1 0.5],'LineWidth',0.005,'EdgeColor','k','LineStyle','--')

plot([0 14.5],[0 0],'k')

for ii=1:4
    s(ii)=plot(1:14,diff(1,:,ii),'-*','lineWidth',4,'Color',c3(ii,:));
end
hold off
set(gca,'FontSize',15)
legend(s,{'Chemophysical','Water','Atrophy','Iron'},'Color',[1 1 1 0.5],'EdgeColor','none','FontSize',18)
set(gca, 'xtick',1:14, 'xticklabel', '','FontSize',15);
ylabel('Effect size','fontSize',18);
xlim([0.5 14.5]);
ax = gca;
ax.YGrid = 'on';
ylim([-1.5 3]);
set(gcf, 'Position',[1 1 1053 618]);


end