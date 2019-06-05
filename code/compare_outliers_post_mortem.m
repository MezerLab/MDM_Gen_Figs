function [mdl]=compare_outliers_post_mortem(color,slopes_score,score,mdm_coef)


% this function generates sup. fig. 10C: The similarity between the molecular variability and the MDM variability in the same brain after removing outliers. 

%% outliers
 % Get the mean and standard deviation of the vector
theMean = mean(score);
stdDev = std(score);

% Get a logical vector of the locations where the value is more than 3 sd away from the mean.
locationsAwayFromMean = abs(score - theMean) > 1.5*stdDev;
outli2=find(locationsAwayFromMean);

% remove outliers
MTV_MRI_out=slopes_score;
MTV_histology_out=score;
MTV_MRI_out(outli2)=[];
MTV_histology_out(outli2)=[];

% compute correlations:
mdl = fitlm(MTV_MRI_out,MTV_histology_out);
slopes_score=MTV_MRI_out;
score=MTV_histology_out;
 
%% plot figure

% find CI:
cifig=figure;
h=plot(mdl);
CI(:,1)=h(4).YData;
CI(:,2)=h(3).YData;
Xd=h(3).XData;
close(cifig)

%figure
% plot fitted line:
linvec=[min(slopes_score),max(slopes_score)];
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
plot(linvec,yprediction,'k--');
hold on
% plot CI:
linePre2=Xd.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
[hl,hp]=boundedline(Xd,linePre2,[linePre2'-CI(:,2),CI(:,1)-linePre2'],'alpha','transparency',0.2,'nan','fill','cmap',color,'orientation','vert') ;
delete(hl);
% plot data points:
for r=1:length(score)
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
% plot signigicance:
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq=['R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f'),' ',stars];
ypos=max(get(gca, 'ylim'));
xpos=min(get(gca, 'xlim'));
text(xpos-0.05*xpos,ypos-0.1*ypos,eq,'FontSize',15)
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