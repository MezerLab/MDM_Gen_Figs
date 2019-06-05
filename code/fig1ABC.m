function fig1ABC(main_dir)

% this fucntion generates Fig.1 A-C and SupFig. 1A-B

%% Bring data

% qMRI parameters:
Par={'R1 [S^-^1]','R2 [S^-^1]','MTsat [p.u.]'};

% phantom lipids type:
lip={'PtdCho','Spg','PS','PtdCho:Chol','PI:PtdCho'};

% load phantom data:
load(fullfile(main_dir,'/fig1/phantom_data1.mat'));

% This .mat file contains the MRI data of the phantoms:
% - LipidMat & LipidMaterror are 55X3X5 matrices contaning mean and std of
%   3 qMRI parameters (see the variable Par) for 5 lipids samples (see the variable lip) with varying water content.
% - MTVMat and MTVMaterror are 55X3X5 matrices contaning mean and std of
%   MTV for 5 lipids samples (see the variable lip) with varying water content.


%% Compute the linear fit between MTV and qMRI parameters for each lipid

% filter MTV values:
MTVMat(MTVMat<0.05)=nan;

% fit the linear fit between MTV and qMRI parameters for each lipid
fitobjectT=[];
slopeMat=[];
for jj=1:length(lip)
    for ii=1:length(Par)
        [v1,I]=sort(MTVMat(:,ii,jj));
        v2=LipidMat(:,ii,jj);
        v2=v2(I);
        v2=v2(~isnan(v1));
        v1=v1(~isnan(v1));
        fitobjectT{ii,jj} = fit(v1,v2,'poly1');
        slopeMat(jj,ii)=fitobjectT{ii,jj}.p1;        
    end
end

%% Plot Fig 1A-B & SupFig. 1A: The dependency of qMRI parameters on the lipid concentration (MTV) for two lipids

% set colors
color_earth=[174  182 185; 122 131  144;...
68  82  93; 174  201  141; 104  136  67;...
51  60  30; 241  219  109;...
192  116  6; 104 50  12; 185  153  112;...
105  73  23; 71  40  22;...
152  143  115; 107  91  77; 53  40  31]/255;
colors=color_earth([2,5,10,7,8],:);

% circled data points:
legendInfo=[];
lip_tmpA{1}=[1,2];
lip_tmpA{2}=[1,5];
lip_tmpA{3}=[3,4];
ind{1}=[32,7,4];
ind{2}=[31,5];
ind{3}=[7,17];
ind{2}=[];
ind{3}=[];

% figures
for ii=1:length(Par)
    lip_tmp=lip_tmpA{ii};
    ind_tmp=ind{ii};
    h=figure;
    hold on
    % create legend:
    for jj=1:length(lip_tmp)
        plot([0,0],[0,0],'LineWidth',2,'color',colors(lip_tmp(jj),:));
        legendInfo{jj}=strcat(lip{lip_tmp(jj)});
    end
    % plot data points:
    for jj=1:length(lip_tmp)
        v1=MTVMat(:,ii,lip_tmp(jj));
        v2=LipidMat(:,ii,lip_tmp(jj));
        plot(v1,v2,'.','color',colors(lip_tmp(jj),:),'MarkerSize',10)
        if ~isempty(ind_tmp)
            plot(v1(ind_tmp(jj)),v2(ind_tmp(jj)),'O','color','k','MarkerSize',10)
            if jj==2
                plot(v1(ind_tmp(jj+1)),v2(ind_tmp(jj+1)),'O','color','k','MarkerSize',10)
            end
        end
    end
    % plot 95% confidence bounds:
    for jj=1:length(lip_tmp)
        [v1,I]=sort(MTVMat(:,ii,lip_tmp(jj)));
        v2=LipidMat(:,ii,lip_tmp(jj));
        v2=v2(I);
        v2=v2(~isnan(v1));
        v1=v1(~isnan(v1));
        fitobject = fitobjectT{ii,lip_tmp(jj)};
        p22 = predint(fitobject,v1,0.95,'functional','on');
        lineV=v1*fitobject.p1+fitobject.p2;
        [hl,hp]=boundedline(v1,lineV,[lineV-p22(:,1),p22(:,2)-lineV],'alpha','transparency',0.2,'nan','fill','cmap',colors(lip_tmp(jj),:),'orientation','vert') ;
    end
    % plot fitted line:
    for jj=1:length(lip_tmp)
        fitobject = fitobjectT{ii,lip_tmp(jj)};
        [v1,I]=sort(MTVMat(:,ii,lip_tmp(jj)));
        v1=v1(~isnan(v1));
        plot([0 0.4],[0 0.4]*fitobject.p1 +fitobject.p2,'--','color',colors(lip_tmp(jj),:),'LineWidth',1);
        plot(v1,v1*fitobject.p1 +fitobject.p2,'color',colors(lip_tmp(jj),:),'LineWidth',2);
    end
    hold off
    legend(legendInfo,'Location','northwest')
    xlabel('MTV [fraction]')
    ylabel(Par{ii})
    grid on
    xlim([0.04 0.35]);
    set(gca,'FontSize',14);
    set(gcf, 'Position',[1 1 1453./3 418]);
end

%% plot Fig. 1A inset: the ambiguity in the biological interpretation of R1. 

legendInfo=[];

% choose data points to present:
lip_tmpA{1}=[1,2];
lip_tmp=lip_tmpA{ii};
ind_tmp=ind{ii};
vecMTV(1)=MTVMat(4,1,2);
vecMTV(2)=MTVMat(7,1,2);
vecMTV(3)=MTVMat(32,1,1);
vecR1(1)=LipidMat(4,1,2);
vecR1(2)=LipidMat(7,1,2);
vecR1(3)=LipidMat(32,1,1);
vecMTVe(1)=MTVMaterror(4,1,2);
vecMTVe(2)=MTVMaterror(7,1,2);
vecMTVe(3)=MTVMaterror(32,1,1);
vecR1e(1)=LipidMaterror(4,1,2);
vecR1e(2)=LipidMaterror(7,1,2);
vecR1e(3)=LipidMaterror(32,1,1);

% create horizontal bar plot:
colors_bar=[colors(2,:) ; colors(2,:); colors(1,:) ];

left =vecR1 ; % Scale 0-1, grows leftwards
right =vecMTV; % Scale 0-35, grows rightwards

% Automatically determine the scaling factor using the data itself
scale = max(right) / max(left);
figure
hold on

% Create the left bar by scaling the magnitude
for ii=1:numel(left)
    b(ii)=herrorbarbar(ii, -left(ii) * scale,-vecR1e(ii) * scale);
    set(b(ii), 'FaceColor', colors_bar(ii,:));
    hold on
end
for ii=1:numel(right)
    b(numel(left)+ii)=herrorbarbar(ii, right(ii),vecMTVe(ii));
    
    set(b(numel(left)+ii), 'FaceColor', colors_bar(ii,:));
    hold on
end
% Now alter the ticks.
yticks = get(gca, 'xtick');

% Get the current labels
labels = get(gca, 'xtickLabel');

if ischar(labels)
    labels = cellstr(labels);
end

% Figure out which ones we need to change
toscale = yticks < 0;

% Replace the text for the ones < 0
labels(toscale) = arrayfun(@(x)sprintf('%0.1e', x), ...
    abs(yticks(toscale) / scale), 'uniformoutput', false);

% Update the tick locations and the labels
set(gca, 'xtick', yticks, 'xticklabel', labels)

yticks = get(gca, 'ytick');
set(gca, 'ytick', [1 2 3], 'yticklabel', {'1','2','3'})

% Now add a different label for each side of the x axis
xmax = max(get(gca, 'xlim'));
label(1) = text(xmax / 2,0.2,  'MTV [fraction]');
label(2) = text(-xmax/ 2, 0.2,'R1 [S^-^1]');

ylim([0.5,3.5]);
set(gca,'FontSize',13);
set(label, 'HorizontalAlignment', 'center', 'FontSize', 14)

%% plot Fig.1C: Unique MRI signatures of  brain lipids

% arrange data for spider plot:
SpiderMat=[];
Total=[slopeMat]';
slopeMatErr=zeros(size(slopeMat));
SpiderMat(:,:,1)=Total+[slopeMatErr ]';
SpiderMat(:,:,2)=Total-[slopeMatErr ]';

% determine axes limits:
Multi=[];
Multi(:,1) = min(min(SpiderMat,[],[3]),[],2);
Multi(:,2) = max(max(SpiderMat,[],[3]),[],2);

% lipids to plot:
lip_tmp=[2,3];

% plot figure:
s=figure('Position', get(0, 'Screensize'));
for ii=1:length(lip_tmp)
    opt_axes=[];
    opt_lines=[];
    opt_area=[];
    opt_area.err        = 'std';
    opt_area.FaceAlpha  = 0;
    opt_area.Color      = colors(lip_tmp(ii),:);
    opt_lines.LineWidth = 3;
    opt_lines.LineStyle = '-';
    opt_lines.Marker    = 'none';
    opt_lines.Labels    = true;
    opt_axes.Phantom    = true;
    opt_axes.Multi   = Multi;
    opt_lines.Color     = colors(lip_tmp(ii),:);
    subplot(1,3,ii)
    opt_axes.Background = 'w';
    opt_axes.Labels     = {'dR1/dMTV','dR2/dMTV','dMTsat/dMTV'};
    polygonplotMDM(SpiderMat(:,lip_tmp(ii),:),opt_axes,opt_lines,opt_area);
    
end
ord=[5:-1:1];
opt_axes=[];
opt_lines=[];
opt_area=[];
opt_area.err        = 'std';
opt_area.FaceAlpha  = 0;
opt_area.Color      = colors(ord,:);
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_lines.Labels    = true;
opt_axes.Phantom    = true;
opt_axes.Multi   = Multi;
opt_lines.Color     = colors(ord,:);
subplot(1,3,3)
opt_axes.Background = 'w';
opt_axes.Labels     = {'dR1/dMTV','dR2/dMTV','dMTsat/dMTV'};
subplot(1,3,3)
polygonplotMDM(SpiderMat(:,ord,:),opt_axes,opt_lines,opt_area);

%% plot SupFig. 1B: The dependency of R1, R2 and MTsat on MTV for 5 lipids 

% no circeled data points in this case:
legendInfo=[];
lip_tmpA{1}=[1:5];
lip_tmpA{2}=[1:5];
lip_tmpA{3}=[1:5];
ind{1}=[];
ind{2}=[];
ind{3}=[7,17];
ind{2}=[];
ind{3}=[];

% plot figure:
h=figure;
for ii=1:length(Par)
    subplot(1,3,ii)
    lip_tmp=lip_tmpA{ii};
    ind_tmp=ind{ii};
    hold on
    % create legend:
    for jj=1:length(lip_tmp)
        plot([0,0],[0,0],'LineWidth',2,'color',colors(lip_tmp(jj),:));
        legendInfo{jj}=strcat(lip{lip_tmp(jj)});
    end
    % plot data points:
    for jj=1:length(lip_tmp)
        v1=MTVMat(:,ii,lip_tmp(jj));
        v2=LipidMat(:,ii,lip_tmp(jj));
        plot(v1,v2,'.','color',colors(lip_tmp(jj),:),'MarkerSize',10)
        if ~isempty(ind_tmp)
            plot(v1(ind_tmp(jj)),v2(ind_tmp(jj)),'O','color','k','MarkerSize',10)
            if jj==2
                plot(v1(ind_tmp(jj+1)),v2(ind_tmp(jj+1)),'O','color','k','MarkerSize',10)
            end
        end
    end
    % plot 95% confidence interval:
    for jj=1:length(lip_tmp)
        [v1,I]=sort(MTVMat(:,ii,lip_tmp(jj)));
        v2=LipidMat(:,ii,lip_tmp(jj));
        v2=v2(I);
        v2=v2(~isnan(v1));
        v1=v1(~isnan(v1));
        fitobject = fitobjectT{ii,lip_tmp(jj)};
        p22 = predint(fitobject,v1,0.95,'functional','on');
        lineV=v1*fitobject.p1+fitobject.p2;
        [hl,hp]=boundedline(v1,lineV,[lineV-p22(:,1),p22(:,2)-lineV],'alpha','transparency',0.2,'nan','fill','cmap',colors(lip_tmp(jj),:),'orientation','vert') ;
    end
    % plot fitted line:
    for jj=1:length(lip_tmp)
        fitobject = fitobjectT{ii,lip_tmp(jj)};
        [v1,I]=sort(MTVMat(:,ii,lip_tmp(jj)));
        v1=v1(~isnan(v1));
        plot([0 0.4],[0 0.4]*fitobject.p1 +fitobject.p2,'--','color',colors(lip_tmp(jj),:),'LineWidth',1);
        plot(v1,v1*fitobject.p1 +fitobject.p2,'color',colors(lip_tmp(jj),:),'LineWidth',2);
    end
    hold off
    xlabel('MTV [fraction]')
    ylabel(Par{ii})
    grid on
    xlim([0.04 0.35]);
    set(gca,'FontSize',14);
end
legend(legendInfo,'Location','northwest')
set(gcf, 'Position',[1 1 1453 418]);

end