function supfig2A(main_dir)

% this function generates SupFig. 2A:  the linearity of the mixtures

%% data

% load lipid mixtures data data:
load(fullfile(main_dir,'/fig1/phantom_data3.mat'));
% this .mat file contains:
% - datay is a cell with 3 vectors each corresponding to a different qMRI
%   parameter (R2, R1, MTsat).
% - MTV is a vector with the MTV values.
%   fitted_line_x, fitted_line_y are cells with 3 vectors each corresponding to a fit between a different qMRI
%   parameter (R2, R1, MTsat) and MTV.
% - R2 is a vector with the R^2 of the fit between different qMRI parameters (R2, R1, MTsat) and MTV.


% set color
colr=[129 187 184]./255 ;
%% plot figure R2

h=figure;
subplot(1,3,1)
plot(fitted_line_x{1},fitted_line_y{1},'--','Color','k','LineWidth',1);
hold on
scatter(MTV, datay{1},100,colr,'filled');
eq={strcat('R^2=',num2str(R2(1),'%.2g'))};
xlabel('MTV [fraction]'); ylabel('R2 [S^-^1]');
set(gca,'FontSize',14);
annotation('textbox',[0.14 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14);
grid on
xticks([0 0.1 0.2 0.3 0.4 0.5]);


%% plot figure  R1
subplot(1,3,3)
plot(fitted_line_x{2},fitted_line_y{2},'--','Color','k','LineWidth',1);
hold on
scatter(MTV, datay{2},100,colr,'filled');
eq={strcat('R^2=',num2str(R2(2),'%.2g'))};
xlabel('MTV [fraction]'); ylabel('R1 [S^-^1]');
set(gca,'FontSize',14);
annotation('textbox',[0.14*5 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14);
grid on
xticks([0 0.1 0.2 0.3 0.4 0.5]);
%% plot figure  MTsat
subplot(1,3,2)
plot(fitted_line_x{3},fitted_line_y{3},'--','Color','k','LineWidth',1);
hold on
scatter(MTV, datay{3},100,colr,'filled');
eq={strcat('R^2=',num2str(R2(3),'%.2g'))};
xlabel('MTV [fraction]'); ylabel('MTsat [p.u.]');
set(gca,'FontSize',14);
annotation('textbox',[0.14*3 0.77 0.3 0.15],'EdgeColor','none','String',eq,'FontSize',14);
grid on
xticks([0 0.1 0.2 0.3 0.4 0.5]);
set(gcf, 'Position',[1 1 1453 418]);

end