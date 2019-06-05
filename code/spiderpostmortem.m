
function spiderpostmortem(MDM,roi,color_new)

% this function plots sup. fig. 10B: Unique MDM signatures of 30 brain regions in two porcine brains. 

% set axes limits:
Multi(:,1)=squeeze(min(MDM(:,1,:)));
Multi(:,2)=squeeze(max(MDM(:,1,:)));

% plot spider plots:
color_spider=color_new;
opt_axes=[];
opt_lines=[];
opt_area=[];
opt_area.FaceAlpha  = 0.5;
opt_area.Color      = color_spider(roi,:);
opt_lines.LineWidth = 3;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_axes.Labels     = {'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'};
opt_lines.Legend    = [];
opt_axes.Multi   = Multi;
opt_lines.Color     = color_spider(roi,:);
opt_axes.Background = 'w';
set(get(gca,'title'),'Position',[-1 1.0060 0])
set(get(gca,'title'),'FontSize',12)
polygonplotMDM((squeeze(MDM(roi,1,:)))',opt_axes,opt_lines,opt_area);
end