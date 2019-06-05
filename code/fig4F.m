
function fig4F(main_dir)

% this function generates fig. 4F: The dependency of MTsat on MTV in three example brain regions. 

%% arragne colors
color_new=[141,211,199;...
255,255,179;...
190,186,218;...
251,128,114;...
128,177,211;...
217,217,217;...
179,222,105;...
188,128,189;...
253,180,98;...
252,205,229;...
204,235,197;...
255,237,111;...
106,61,154;...
51,160,44;...
31,120,180]./255;

%% Example ROIs
thal=[153 120 175]./255;
wm=  [ 5 149 198]./255;
pons=[34 100 95]./255;

roi=[7,8,9];
color_new(7,:)=wm;
color_new(8,:)=pons;
color_new(9,:)=thal;
ROI_list={'WM','Pons','Thalamus'};

%% load MRI data 

x=load(fullfile(main_dir,'/fig4/porcine_brain2_additional_data.mat'));

Parameter_range=[2.3 5.6]; % range of MTsat values
TV_range=[0.19 0.37]; % range of TV values
Parameter_str='MTsat [p.u.]';

BinTVm=x.BinTVm; % median values of MTV bins used for the linear fit
BinParm=x.BinParm; % median values of MTsat bins used for the linear fit
STDm=x.STDm; % std in each MTsat bin
fit=x.fit; % fit information

%% plot fig. 4F

TV_str='MTV [fraction]';

hold on
% prepare legend
for i=1:length(roi)
    plot([0,0],[0,0],'LineWidth',2,'color',color_new(roi(i),:));
    legendInfo{i}=strcat(ROI_list{(i)});
end
% plot data points
for i=1:length(roi)
    v1=BinTVm{roi(i)};
    v2=BinParm{roi(i)};
    plot(v1,v2,'.','color',color_new(roi(i),:),'MarkerSize',15)
end
% plot std in each bin
for i=1:length(roi)
    v1=BinTVm{roi(i)};
    err=STDm{roi(i)};
    [hl,hp]=boundedline(v1,v1*fit(roi(i),1)+fit(roi(i),2),err,'alpha','transparency',0.2,'cmap',color_new(roi(i),:)) ;
end
% plot fitted line
for i=1:length(roi)
    v1=BinTVm{roi(i)};
    plot(TV_range,TV_range*fit(roi(i),1) +fit(roi(i),2),'--','color',color_new(roi(i),:),'LineWidth',1);
    plot(v1,v1*fit(roi(i),1) +fit(roi(i),2),'color',color_new(roi(i),:),'LineWidth',2);
end
hold off
legend(legendInfo)
xlim(TV_range)
ylim(Parameter_range);
xlabel(TV_str)
ylabel(Parameter_str)
set(gca,'FontSize',15);
grid on
set(gca,'Color',[0.93 0.93 0.93])

end

