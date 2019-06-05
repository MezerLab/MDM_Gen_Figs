
function fig2d_spider(color_new,comrl,young_ind,fit_info,v,l,r)

% this function generates figure 2D: Unique MDM signatures for different brain regions

%% arrange data

% group young subjects' MDM measurements for left ROIs:
Test_young=nan(4,length(v),length(young_ind).*2);
Test_young(:,:,1:length(young_ind))=[rot90(fit_info.data{1}(v,[1],young_ind)) ; rot90(fit_info.data{2}(v,[1],young_ind)); rot90(fit_info.data{3}(v,[1],young_ind)); rot90(fit_info.data{4}(v,[1],young_ind))];

% joit right ROIs as well:
for ii=1:length(v)
    if ismember(v(ii),l) && comrl
        ind=find(l==v(ii));
        Test_young(:,ii,length(young_ind)+1:end)=[rot90(fit_info.data{1}(r(ind),[1],young_ind)) ; rot90(fit_info.data{2}(r(ind),[1],young_ind)); rot90(fit_info.data{3}(r(ind),[1],young_ind)); rot90(fit_info.data{4}(r(ind),[1],young_ind))];
    end
end

%% Figure 2D

%  spiderplots: 
Multi=[];
color_spider=color_new;

ROI_list_new={'White-Matter','Thalamus','Putamen','Pallidum','Hippocampus','Frontal Cortex'};

% set axes limits:
Multi(:,1) = min(nanmedian(Test_young,3),[],[2]);
Multi(:,1) = Multi(:,1)-0.5*abs(Multi(:,1));
Multi(:,2) = max(nanmedian(Test_young,3),[],[2]);
Multi(:,2) = Multi(:,2)+0.1*abs(Multi(:,2));

% plot figure
xpos=0.13:0.09:0.13*5;
ypos=0.8007:-0.1457:0.8007-0.1457*2;
s=figure('Position', get(0, 'Screensize'));
for ii=1:length(v)
    opt_axes=[];
    opt_lines=[];
    opt_area=[];
    opt_area.err        = 'std';
    opt_area.FaceAlpha  = 0.5;
    opt_area.Color      = color_spider{ii};
    opt_lines.LineWidth = 2;
    opt_lines.LineStyle = '-';
    opt_lines.Marker    = 'none';
    opt_lines.Labels    = false;
    opt_lines.Legend    = [];
    opt_axes.Multi   = Multi;
    opt_lines.Color     = color_spider{ii};
    opt_axes.Background = 'w';
    opt_axes.subPlot=true;
    subplot(5,3,ii)
    set(get(gca,'title'),'Position',[-1 1.0060 0])
    set(get(gca,'title'),'FontSize',12)
    polygonplotMDM((Test_young(:,ii,:)),opt_axes,opt_lines,opt_area);
    pos = get(gca, 'Position');
    [x,y] = ind2sub([5,3],ii);
    pos(1)=xpos(x);
    pos(2)=ypos(y);
    set(gca, 'Position', pos)
end
color_spider=color_new;
color_spider=cell2mat(color_spider);
ord=[7:14,1:6];
opt_axes=[];
opt_lines=[];
opt_area=[];
opt_area.err        = 'std';
opt_area.FaceAlpha  = 0.5;
opt_area.Color      = color_spider(ord,:);
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_lines.Labels    = false;
opt_lines.Legend    = [];
opt_axes.Multi   = Multi;
opt_lines.Color     = color_spider(ord,:);
opt_axes.Background = 'w';
opt_axes.subPlot=true;
subplot(5,3,15)
set(get(gca,'title'),'Position',[-1 1.0060 0])
set(get(gca,'title'),'FontSize',12)
polygonplotMDM((Test_young(:,ord,:)),opt_axes,opt_lines,opt_area);
pos = get(gca, 'Position');
[x,y] = ind2sub([5,3],ii);
pos(1)=xpos(5);
pos(2)=ypos(3);
set(gca, 'Position', pos)
set(gcf, 'Position',[1 1 1140 500]);

%% Axis legend
figure;
opt_axes=[];
opt_lines=[];
opt_area=[];
opt_area.err        = 'std';
opt_area.FaceAlpha  = 0.5;
opt_area.Color      = [ 1 1 1; 1 1 1 ];
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_lines.Legend    = [];
opt_lines.Labels    = true;
opt_axes.Multi   = Multi;
opt_lines.Color     = [ 1 1 1; 1 1 1 ];
opt_axes.Background = 'w';
opt_axes.Labels={'dMTsat','dR1','dMD','dR2'};
tmp=[];
tmp(:,:,1)=Multi;
tmp(:,:,2)=Multi;
set(get(gca,'title'),'FontSize',18)
polygonplotMDM((tmp),opt_axes,opt_lines,opt_area);
set(gcf, 'Position',[1 1 693 208]);
end