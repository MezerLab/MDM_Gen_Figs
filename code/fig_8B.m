
function fig_8B(diff_M)

% this function generates figure 7B: The spatial correlation between aging-related changes in different biological markers. 

% color:
c=[156 163 169; 161 162 120 ; 194 172 124; 152 136 146]./255;


% plot figure 7B

lables={'Chemophysical Effect size','Water Effect size','Volume Effect size','Iron Effect size'};

% order of comparisons between markers:
xpar=[1,1,2,2,1,3];
ypar=[3,4,3,4,2,4];

pvals=[];
mdlL=[];

% compute correlations between all aging markers:
for ii=1:length(xpar)
    x=diff_M(:,xpar(ii));
    y=diff_M(:,ypar(ii));
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
    linvec=[min(x) max(x)];
    linePre=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1}; % fitted line
    linePre2=Xd.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1}; % prediction for each data point
    
    % save fit info:
    boundL{ii}=[linePre2'-CI(:,2),CI(:,1)-linePre2']; % CI for each data point
    XdL{ii}=Xd;
    linePre2L{ii}=linePre2; 
    linePreL{ii}=linePre;
    linvecL{ii}=linvec;
    mdlL{ii}=mdl;
    pvals(ii)=mdlL{ii}.Coefficients{2,4};
end


% plot figure:
figure(2);
for ii=1:length(xpar)
    subplot(3,2,ii)
    hold on
    x=diff_M(:,xpar(ii));
    y=diff_M(:,ypar(ii));
    [hl,hp]=boundedline(XdL{ii},linePre2L{ii},boundL{ii},'alpha','transparency',0.2,'nan','fill','cmap',[0.5 0.5 0.5],'orientation','vert') ;
    delete(hl);
    scatter(x,y,20,[0.5 0.5 0.5],'filled');
    plot(linvecL{ii},linePreL{ii},'--', 'LineWidth',2,'Color',[0.5 0.5 0.5])
    [stars]=pval2stars(pvals(ii),'num');
    eq=['R^2=',num2str(mdlL{ii}.Rsquared.Adjusted,'%0.2g'),', ', stars];
    xlabel(lables{xpar(ii)});
    ylabel(lables{ypar(ii)});
    title(eq)
    grid on
    hold off
    ax = gca;
    ax.XColor=c(xpar(ii),:);
    ax.YColor=c(ypar(ii),:);
    ax.FontWeight='bold';
end
set(gcf, 'Position',[1 1 1053 800]);

end