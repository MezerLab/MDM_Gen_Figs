function fig1D(main_dir)

% This function generates Fig. 1D and SupFig. 1E-G
%% Bring data

% load phantom data
load(fullfile(main_dir,'/fig1/phantom_data2.mat'));

% This .mat file contains:
% - mix1,2,3 are 2X7 matrices with the rows representing predicted and measured MTsat values and the columns are different water concentrations for 3 mixtures with
%   different PtdCho to PS ratios.
% - xxR1,xxR2,xxMTsat and yyR1,yyR2,yyMTsat are the predicted and real
%   R1, R2 and MTsat values for all lipid samples.

%% set colors

color_earth=[174  182 185; 122 131  144;...
68  82  93; 174  201  141; 104  136  67;...
51  60  30; 241  219  109;...
192  116  6; 104 50  12; 185  153  112;...
105  73  23; 71  40  22;...
152  143  115; 107  91  77; 53  40  31]/255;


colr=[129 187 184]./255 ;

%% Fig. 1D: Predicting the MRI signal of a lipid mixture from the signal of pure lipids. 

h=figure;
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle='--';
hold on
m(1)=scatter(mix1(1,:),mix1(2,:),100,color_earth(13,:),'filled','LineWidth',1.2);
m(2)=scatter(mix2(1,:),mix2(2,:),100,color_earth(14,:),'filled','LineWidth',1.2);
m(3)=scatter(mix3(1,:),mix3(2,:),100,color_earth(15,:),'filled','LineWidth',1.2);

xx=[mix1(1,:) mix2(1,:) mix3(1,:)];
yy=[mix1(2,:) mix2(2,:) mix3(2,:)];

mdl = fitlm(xx,yy);
R2=mdl.Rsquared.Adjusted;
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(R2,'%.2g')),stars};
xlabel('MTsat [p.u.]'); ylabel('MTsat predicted [p.u.]')
set(gca,'FontSize',15)
l=legend(m,'PC:PS(1:1)','PC:PS(2:1)','PC:PS(1:2)','Location','southeast');
l.FontSize=14;
annotation('textbox',[0.14 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14)
xticks([0 0.2 0.4 0.6 0.8]);
yticks([0 0.2 0.4 0.6 0.8]);
grid on
xlim([0.07 0.7]); ylim([0.07 0.7]);

%% SupFig. 1C: Prediction for the R1, R2 and MTsat of 9 lipid mixtures (y-axes) vs. their true values (x-axes)

h=figure;
subplot(1,3,3);
scatter(xxMT,yyMT,60,colr,'filled');
hold on
identityLine(gca);
mdl = fitlm(xxMT,yyMT);
R2=mdl.Rsquared.Adjusted;
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(R2,'%.2g')),stars};
xlabel('MTsat [p.u.]'); ylabel('MTsat predicted [p.u.]');
set(gca,'FontSize',14);
annotation('textbox',[0.14*5 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14);
xlim([0 1.65]); ylim([0 1.65]);

subplot(1,3,1);
scatter(xxR1,yyR1,60,colr,'filled');
hold on
mdl = fitlm(xxR1,yyR1);
R2=mdl.Rsquared.Adjusted;
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(R2,'%.2g')),stars};
identityLine(gca);
xlabel('R1 [S^-^1]'); ylabel('R1 predicted [S^-^1]');
set(gca,'FontSize',14);
annotation('textbox',[0.14 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14);

subplot(1,3,2);
scatter(xxR2,yyR2,60,colr,'filled');hold on
hold on
mdl = fitlm(xxR2,yyR2);
R2=mdl.Rsquared.Adjusted;
[stars]=pval2stars(mdl.Coefficients{2,4},'num');
eq={strcat('R^2=',num2str(R2,'%.2g')),stars};
identityLine(gca);
xlabel('R2 [S^-^1]'); ylabel('R2 predicted [S^-^1]');
set(gca,'FontSize',14);
annotation('textbox',[0.14*3 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14);

set(gcf, 'Position',[1 1 1453 418]);

end