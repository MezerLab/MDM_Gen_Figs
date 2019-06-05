function [mdl,slopes_score,score]=compare_lipids_post_mortem(ROIs_labels,MDM,lipids,lipids_labels,mdm_coef,plot_indiv_lipids,color)

% This function generates figure 4H,E,D: The MDM signatures explain the molecular variability across brain regions. Standard qMRI parameters do not explain the molecular variability

% Compute 1st PC of lipidomics:

[coefflipids,scorelipids,latent,tsquared,explainedlipids,mu] = pca(lipids);
score=scorelipids(:,1);

% Compute 1st PC of MRI:
% for mdm_coef=1: Compute 1st PC of MDM:
% for mdm_coef=4: Compute 1st PC of standard qMRI parameters:

slopes=squeeze(MDM(:,mdm_coef,[1,2,3,4]));
[coeff,slopes_score,latent,tsquared,explained,mu] = pca(zscore(slopes));
slopes_score=slopes_score(:,1);

% for mdm_coef=3: use MTV data instead:
if mdm_coef==3
    slopes_score=MDM(:,mdm_coef,2);
end

% compute correlation between MRI and histology:
mdl = fitlm(slopes_score,score);

% find CI:
cifig=figure;
h=plot(mdl);
CI(:,1)=h(4).YData;
CI(:,2)=h(3).YData;
Xd=h(3).XData;
close(cifig)

%% plot figure

% plot fitted line
linvec=[min(slopes_score),max(slopes_score)];
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
plot(linvec,yprediction,'k--');

% plot CI:
hold on
linePre2=Xd.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
[hl,hp]=boundedline(Xd,linePre2,[linePre2'-CI(:,2),CI(:,1)-linePre2'],'alpha','transparency',0.2,'nan','fill','cmap',color,'orientation','vert') ;
delete(hl);
% plot data points:
for r=1:length(ROIs_labels)
    h(r)=scatter(slopes_score(r),score(r),100,color,'filled');
end
grid on
hold off

if mdm_coef==1
    xlabel('1st PC MDM');
elseif mdm_coef==4
    xlabel('1st PC standard qMRI');
elseif mdm_coef==3
    xlabel('MTV [fraction]');
end

ylabel('1st PC Lipidomics');
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
% plot significance:
eq=['R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f'),' ',stars];
ypos=max(get(gca, 'ylim'));
xpos=min(get(gca, 'xlim'));
if mdm_coef==1
text(xpos-0.9*xpos,ypos-0.7*ypos,eq,'FontSize',15)
else
text(xpos-0.1*xpos,ypos-0.1*ypos,eq,'FontSize',15)
end
set(gca,'FontSize',15)
NumTicks = 6;
L = get(gca,'YLim');
xlim([-2.55 4.25]);
if mdm_coef==3
    xlim([0.12 0.35]);
    
end
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gcf, 'Position',[1 1 1268 588]);
set(gca,'Color',[0.93 0.93 0.93])

end