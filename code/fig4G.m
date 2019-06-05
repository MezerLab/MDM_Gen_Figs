
function fig4G(MDM)
% this function plots figure 4G: Unique MDM signatures for different brain regions. 

% remove data of the other porcine brain:
MDM=MDM(16:end,:,:);

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


%% plot fig. 4G

TV_str='MTV [fraction]';

% set axes limits:
Multi(:,1)=squeeze(min(MDM(:,1,:)));
Multi(:,2)=squeeze(max(MDM(:,1,:)));

% spider plot:
color_spider=color_new;
opt_axes=[];
opt_lines=[];
opt_area=[];
opt_area.FaceAlpha  = 0.5;
opt_area.Color      = color_spider(roi,:);
opt_lines.LineWidth = 4;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_axes.Labels     = {'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'};
opt_lines.Legend    = [];
opt_axes.Multi   = Multi;
opt_lines.Color     = color_spider(roi,:);
opt_axes.Background = 'w';
polygonplotMDM((squeeze(MDM(roi,1,:)))',opt_axes,opt_lines,opt_area);
set(gca,'FontSize',15);

end
