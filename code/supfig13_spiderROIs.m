
function supfig13_spiderROIs(color_spider,FDRval,fit_info,Test_young,Test_old,v)

% this function generates SupFig. 13: aging-related changes in MDM signatures
xpos=0:0.16:0+5*0.16;
ypos=0.8007:-0.21:0.8007-0.21*2;
s=figure('Position', get(0, 'Screensize'));
for ii=1:length(v) % loop over ROIs
    subplot(5,3,ii)
    curr=[];
    opt_axes=[];
    opt_lines=[];
    opt_area=[];
    curr(:,:,1)=squeeze(Test_young(:,ii,:));
    curr(:,:,2)=squeeze(Test_old(:,ii,:));
    margins= [curr(:,:,1)  curr(:,:,2)];
    curr=permute(curr,[1 3 2]);
    % set axes limits:
    Multi(:,1)= prctile(margins',5)';
    Multi(:,2)= prctile(margins',95)';
    opt_area.err        = 'std';
    opt_area.FaceAlpha  = 0.5;
    opt_area.Color      = [color_spider(ii,:); [175,175,175]./255];
    opt_lines.Color     = [color_spider(ii,:); [175,175,175]./255];
    opt_lines.LineWidth = 2;
    opt_lines.LineStyle = '-';
    opt_lines.Marker    = 'none';
    opt_lines.Labels    = false;
    opt_axes.subPlot=true;
    opt_axes.Multi= Multi;
    opt_axes.Background = 'w';
    opt_axes.Signif=squeeze(FDRval(1,ii,[1,2,3,4]));
    opt_axes.Labels     = {[fit_info.str{1} ' s'],[fit_info.str{2} ' s'],[fit_info.str{3} ' s'],[fit_info.str{4} ' s']};
    polygonplotMDM(curr,opt_axes,opt_lines,opt_area);
    pos = get(gca, 'Position');
    [x,y] = ind2sub([5,3],ii);
    pos(1)=xpos(x);
    pos(2)=ypos(y);
    pos([3,4])=1.3*pos([3,4]);
   set(gca, 'Position', pos)
end

% show axes labels
subplot(5,3,15)
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
opt_lines.Labels    = false;
opt_lines.Color     = [ 1 1 1; 1 1 1 ];
opt_axes.Background = 'w';
opt_axes.Labels={'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'};
tmp=[];
tmp(:,:,1)=Multi;
tmp(:,:,2)=Multi;
polygonplotMDM((tmp),opt_axes,opt_lines,opt_area);
pos = get(gca, 'Position');
[x,y] = ind2sub([5,3],15);
pos(1)=xpos(5);
pos(2)=ypos(3);
pos([3,4])=1.3*pos([3,4]);
set(gca, 'Position', pos)
end