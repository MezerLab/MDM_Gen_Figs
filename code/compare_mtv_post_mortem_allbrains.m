function [mdl]=compare_mtv_post_mortem_allbrains(MTV_histology,MTV_MRI,color)

% this function generates figure 4C: Validation of the qMRI estimation of MTV.  


%% compute correlation of MRI MTV and histology MTV:

mdl = fitlm(MTV_MRI,MTV_histology);

% find CI:
cifig=figure;
h=plot(mdl);
CI(:,1)=h(4).YData;
CI(:,2)=h(3).YData;
Xd=h(3).XData;
close(cifig)

%% plot figure 

% plot fitted line:
linvec=[min(MTV_MRI),max(MTV_MRI)];
yprediction=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
plot(linvec,yprediction,'k--');

% plot CI:
hold on
linePre2=Xd.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1};
[hl,hp]=boundedline(Xd,linePre2,[linePre2'-CI(:,2),CI(:,1)-linePre2'],'alpha','transparency',0.2,'nan','fill','cmap',color,'orientation','vert') ;
delete(hl);

% plot data points:
for r=1:length(MTV_MRI)
    h(r)=scatter(MTV_MRI(r),MTV_histology(r),100,color,'filled');
end
grid on
hold off
xlabel('MTV MRI');
ylabel('MTV Histology');

% plot significance:
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(mdl.Rsquared.Adjusted,'%0.2f')),[stars]};
ypos=max(get(gca, 'ylim'));
xpos=min(get(gca, 'xlim'));
text(xpos-0.1*xpos,ypos-0.2*ypos,eq,'FontSize',15)
set(gca,'FontSize',15) 
ylim([0.12 0.33]);
set(gca,'Color',[0.93 0.93 0.93])

end
