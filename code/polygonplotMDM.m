
function polygonplotMDM(data,opt_axes,opt_lines,opt_area)

% This is an adaptation of the function polygonplot by Vctor Martnez
% Cagigal.
% ----------------------------------------------------------------------- %
% Function polygonplot plots a kind of radial plot whose shape is a N-poly%
% gon, depending on the size of the data. It have also provide users the  %
% functionality of plotting a shaded error area if the dimension of the   %
% data is highly enough.                                                  %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix.                                          %
%           * Mean case: [N x nlines], with N corresponding to the number %
%             of points of each line (number of polygon sides) and nlines %
%             corresponding to the number of lines to plot.               %
%           * Error bar case: [N x nlines x M], with M corresponding to   %
%             the number of observations.                                 %
%       - opt_axes: (Optional) Struct with the axes options.              %
%           * opt_axes.Ticks:       Vector that contains the axis ticks.  %
%           * opt_axes.Background:  Background color of the tick labels.  %
%           * opt_axes.Labels:      {1 x N} cell-vector with the axes lbls%
%       - opt_lines:(Optional) Struct with the lines options.             %
%           * opt_lines.LineWidth:  Width of the means.                   %
%           * opt_lines.LineStyle:  Line style of the means.              %
%           * opt_lines.Marker:     Marker of the means.                  %
%           * opt_lines.Ccolor:      [nlines x 3] matrix with RGB colors.  %
%           * opt_lines.Labels:     Boolean. If true, text labels of each %
%                                   point are plotted.                    %
%           * opt_lines.Legend:     {1 x nlines} cell-matrix with the     %
%                                   legend of each line.                  %
%       - opt_area: (Optional) Struct with the shaded area options.       %
%           * opt_area.err:         Type of error to plot:                %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
%           * opt_area.FaceAlpha:   Alpha transparency constant.          %
%           * opt_area.Color:       [nlines x 3] matrix with RGB colors.  %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       d_ex = [2 3; 1 0; 0.1 3; -1 7; -0.2 0.9];                         %
%       data = cat(3,d_ex-0.5,d_ex,d_ex+0.7);                             %
%       polygonplot(data);                                                %
% ----------------------------------------------------------------------- %
%   Author:  Vctor Martnez Cagigal                                        %
%   Date:    21/03/2017                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %


% Defaults
if(nargin<4)
    opt_area.err = 'std';
    opt_area.FaceAlpha = 0.5;
end
if(nargin<3)
    opt_lines.LineWidth = 2;
    opt_lines.LineStyle = '-';
    opt_lines.Marker    = 'none';
end
if(nargin<2)
    opt_axes = [];
    opt_axes.Background = 'none';
end
if(nargin<1)
    error('Not enough parameters');
end

% Error detection
if(length(size(data))>3)
    error('Data must be a MxNxnlines or Nxnlines matrix.');
elseif(length(size(data))==3)       % Shaded
    [N,nlines,M] = size(data);
elseif(length(size(data))==2)       % Only mean
    [N,nlines]   = size(data);
    M = 0;
end
if(~isfield(opt_lines,'Color'))     % Color properties
    line_color = [];
else
    if(size(opt_lines.Color,1)~=nlines)
        error('Number of colors must be equal to the number of lines to plot.');
    else
        if(~isfield(opt_area,'Color') && M~=0)
            error('Please, specify also the colors of the shaded areas.');
        else
            if(size(opt_area.Color,1)~=nlines)
                error('Number of shaded areas must be equal to the number of lines to plot.');
            else
                line_color = opt_lines.Color;
                opt_lines = rmfield(opt_lines,'Color');
            end
        end
    end
end
if(~isfield(opt_lines,'Labels'))    % Text labels
    labels = false;
else
    labels = opt_lines.Labels;
    opt_lines = rmfield(opt_lines,'Labels');
end
if(~isfield(opt_axes,'Phantom'))    % Text labels
    Phantom = false;
else
    Phantom = true;
end
if(~isfield(opt_axes,'Multi'))    % Text labels
    Multi = false;
else
    Multi = opt_axes.Multi;
end
if(~isfield(opt_axes,'subPlot'))    % Text labels
    subPlot = false;
else
    subPlot = true;
end
if(~isfield(opt_area,'Color'))
    opt_area.Color = 0.5.*ones(nlines,3);
end
if(~isfield(opt_area,'Confint'))
    Confint = false;
else
    Confint = true;
end
if(~isfield(opt_axes,'Signif'))
    pval=0;
else
    pval=1;
    pvalm=opt_axes.Signif;
end
if(~isfield(opt_area,'FaceAlpha'))
    opt_area.FaceAlpha = 0.5;
end
if(~isfield(opt_lines,'Legend'))
    leg = [];
else
    leg = opt_lines.Legend;
    opt_lines = rmfield(opt_lines,'Legend');
end
if(isfield(opt_axes,'Labels'))
    if(length(opt_axes.Labels)~=N)
        error('You must provide N axis labels.');
    end
end
set(gcf,'Color',[1 1 1]);

% Compute the isocurve Ticks
data_mad  = mad(data,1,3);
data_median = nanmedian(data,3);
totp=data_median+0.5.*data_mad;
totn=data_median-0.5.*data_mad;

%data_min=min(totn,[],2);

data_min = min(min(data,[],[3]),[],2);
% data_max=max(totp,[],2);

data_max = max(max(data,[],[3]),[],2);

if Multi
    data_min=Multi(:,1);
    data_max=Multi(:,2);
end

if Phantom
    data_minN=data_min-0.1*(data_max-data_min);
    % data_maxN=data_max+0.1*(data_max-data_min);
    data_min=data_minN;
    % data_max=data_maxN;
end
if(~isfield(opt_axes,'Ticks'))
    for ii=1:length(data_max)
        opt_axes.Ticks{ii} = linspace(data_min(ii),data_max(ii),9);
    end
elseif(length(opt_axes.Ticks)<2)
    error('There should be more than 2 isocurves.');
end
for ii=1:length(data_max)
    % r_iso(:,ii) = opt_axes.Ticks{ii}(:);
    
    r_iso(:,ii) = (opt_axes.Ticks{ii}(:)-data_min(ii))/(data_max(ii)-data_min(ii));
end
th_iso = (2*pi/N)*(ones(length(opt_axes.Ticks{ii}),1)*(N:-1:1));
[x,y] = pol2cart(th_iso, r_iso);
h_iso = line([x,x(:,1)]',[y,y(:,1)]','LineWidth',0.5,'Color',0.85.*ones(1,3));
uistack(h_iso,'bottom');

for iso_id = 1:1:N   % Exclude axes from legend
    set(get(get(h_iso(iso_id),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold on;
uistack(h_iso,'bottom');
% Compute and plot the axes depending on N
th_jump = (2*pi/N)*(ones(2,1)*(N:-1:1));
radii   = [zeros(1,N); ones(1,N)];
[x,y]   = pol2cart(th_jump, radii);
h_axes  = line(x,y,'LineWidth',1,'Color','k');
for ax_id = 1:1:N   % Exclude axes from legend
    set(get(get(h_axes(ax_id),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold on;
uistack(h_axes,'bottom');

% Display axis Ticks
if(~isfield(opt_axes,'Background')), opt_axes.Background = 'none'; end
for ii=1:length(data_max)
    loc = (opt_axes.Ticks{ii}(:)-data_min(ii))/(data_max(ii)-data_min(ii));
    %  loc = (opt_axes.Ticks(:)-data_min)/(data_max-data_min);
    odd_v = 0.03.*ones(size(opt_axes.Ticks{ii}));
    odd_v(2:2:end) = -0.7.*odd_v(2:2:end);
    %     if ii==1
    %      for t = 1:1:length(opt_axes.Ticks{ii})
    %          text(loc(t),odd_v(t),num2str(opt_axes.Ticks{ii}(t)),'Background',opt_axes.Background,'FontSize',7);
    %      end
    %     end
    %     if ii==2
    %         for t = 1:1:length(opt_axes.Ticks{ii})
    %          text(odd_v(t),-loc(t),num2str(opt_axes.Ticks{ii}(t)),'Background',opt_axes.Background,'FontSize',7);
    %         end
    %     end
    %        if ii==3
    %         for t = 1:1:length(opt_axes.Ticks{ii})
    %          text(-loc(t),odd_v(t),num2str(opt_axes.Ticks{ii}(t)),'Background',opt_axes.Background,'FontSize',7);
    %         end
    %        end
    %           if ii==4
    %         for t = 1:1:length(opt_axes.Ticks{ii})
    %          text(odd_v(t),loc(t),num2str(opt_axes.Ticks{ii}(t)),'Background',opt_axes.Background,'FontSize',7);
    %         end
    %           end
end



% Plot the data
opt_lines_str = adapt_options(opt_lines);
if(M == 0)    % Only mean
    R  = ([data; data(1,:)]-[data_min ; data_min(1,:)])./([data_max; data_max(1)]-[data_min ; data_min(1,:)]);
    TH = (2*pi/N)*((N:-1:0)'*ones(1,nlines));
    [X,Y] = pol2cart(TH, R);
    if(isempty(line_color)), plot(X,Y,opt_lines_str{:});
    else
        for n = 1:1:nlines,  plot(X(:,n),Y(:,n),opt_lines_str{:},'Color',line_color(n,:)); end
    end
    if(labels)
        for n = 1:1:nlines
            for t = 1:1:N
                text(X(t,n),Y(t,n),num2str(data(t,n),'%.1g'),'FontSize',12);
            end
        end
    end
    axis square;
    axis off;
else          % Shaded area
    % Computing the mean and standard deviation of the data matrix
    data_mean = nanmedian(data,3);
    data_std  = mad(data,1,3);
    if Confint
        data_mean=data(:,:,1);
        data_std=data(:,:,[2,3]);
    end
    % Type of error plot
    switch(opt_area.err)
        case 'std', err = data_std;
        case 'sem', err = (data_std./sqrt(size(data,1)));
        case 'var', err = (data_std.^2);
        case 'c95', err = (data_std./sqrt(size(data,1))).*1.96;
    end
    
    % Plots
    if Confint
        m_down = data_mean-err(:,:,1)./2;
        m_up   = data_mean+err(:,:,2)./2;
    else
        m_down = data_mean-err./2;
        m_up   = data_mean+err./2;
    end
    m_down = (m_down-data_min*ones(1,nlines))./(data_max*ones(1,nlines)-data_min*ones(1,nlines));
    m_down = [m_down; m_down(1,:)];
    m_up   = (m_up-data_min*ones(1,nlines))./(data_max*ones(1,nlines)-data_min*ones(1,nlines));
    m_up   = [m_up; m_up(1,:)];
    for ii=1:length(data_min)
        R(ii,:)  = (data_mean(ii,:)-data_min(ii))./(data_max(ii)-data_min(ii));
    end
    R(N+1,:)= (data_mean(1,:)-data_min(1))./(data_max(1)-data_min(1));
    TH = (2*pi/N)*((N:-1:0)'*ones(1,nlines));
    [X,Y] = pol2cart(TH, R);
    
    for n = 1:1:nlines
        [xa,ya] = pol2cart([TH(:,n); fliplr(TH(:,n))],[m_down(:,n); fliplr(m_up(:,n))]);
        pat = fill(xa,ya,opt_area.Color(n,:));
        hold on;
        set(pat, 'EdgeColor', 'none');
        set(pat, 'FaceAlpha', opt_area.FaceAlpha);
        uistack(pat,'top');
        
    end
    if(isempty(line_color))
        h_leg = plot(X,Y,opt_lines_str{:}); hold on;
    else
        for n = 1:1:nlines
            h_leg(n) = plot(X(:,n),Y(:,n),opt_lines_str{:},'Color',line_color(n,:));
        end
    end
    if(labels)
        for n = 1:1:nlines
            for t = 1:1:N
                text(X(t,n),Y(t,n),num2str(data_mean(t,n),'%.3g'),'FontSize',9);
            end
        end
    end
    axis equal;
    axis off;
    if(~isempty(leg))
        leg1=legend(h_leg([1:n]),leg,'Location','northeastoutside');
        set(leg1,'Box','off')
    end
end
% Display axis labels
if(isfield(opt_axes,'Labels'))
    th_lbl = (2*pi/N)*(N:-1:1);
    r_lbl  = 1.1.*ones(1,N);
    [xlbl,ylbl] = pol2cart(th_lbl,r_lbl);
    for lid = 1:1:N
        if subPlot
            if lid==3
                if pval==1 && pvalm(lid)<0.0500
                    [stars]=pval2stars(pvalm(lid),'stars');
                    %   tt(1)=text(xlbl(lid)+0.2*xlbl(lid),ylbl(lid),[opt_axes.Labels{lid}],'FontSize',12);
                    % tt(2)=text(xlbl(lid)+0.2*xlbl(lid),ylbl(lid)-0.1,['p=' num2str(pvalm(lid),'%10.1e')],'FontSize',12);
                    text(xlbl(lid)+xlbl(lid),ylbl(lid),stars,'FontSize',12);
                    
                else
                    %    tt(3)=text(xlbl(lid)+0.2*xlbl(lid),ylbl(lid),opt_axes.Labels{lid},'FontSize',12);
                end
            else
                if pval==1 && pvalm(lid)<0.0500
                    [stars]=pval2stars(pvalm(lid),'stars');
                    %   tt(1)=text(xlbl(lid)+0.2*xlbl(lid),ylbl(lid),[opt_axes.Labels{lid}],'FontSize',12);
                    % tt(2)=text(xlbl(lid)+0.2*xlbl(lid),ylbl(lid)-0.1,['p=' num2str(pvalm(lid),'%10.1e')],'FontSize',12);
                    %  tt(4)=text(xlbl(lid)+0.2*xlbl(lid),ylbl(lid)-0.1,stars,'FontSize',12);
                    text(xlbl(lid),ylbl(lid),[stars],'FontSize',12);
                    
                    %tt(4)=text(xlbl(lid),ylbl(lid),[opt_axes.Labels{lid}],'FontSize',12);
                    %tt(5)=text(xlbl(lid),ylbl(lid)-0.1,['p=' num2str(pvalm(lid),'%10.1e')],'FontSize',12);
                else
                    %   tt(6)=text(xlbl(lid),ylbl(lid),opt_axes.Labels{lid},'FontSize',12);
                end
            end
        elseif Phantom
            if pval==1 && pvalm(lid)<0.0500
                [stars]=pval2stars(pvalm(lid),'stars');
                if lid==4
                    text(xlbl(lid),ylbl(lid),[stars newline opt_axes.Labels{lid} ],'FontSize',12,'HorizontalAlignment','center');
                else
                    text(xlbl(lid),ylbl(lid),[opt_axes.Labels{lid} newline stars],'FontSize',12,'HorizontalAlignment','center');
                end
                
            else
                text(xlbl(lid),ylbl(lid),opt_axes.Labels{lid},'FontSize',12,'HorizontalAlignment','center');
            end
        else
            
            Pos_in=[1.321 -0.045 0; 0 -1.116 0 ;-1.247 -0.045 0 ; 0 1.116 0];
            if pval==1 && pvalm(lid)<0.0500
                [stars]=pval2stars(pvalm(lid),'stars');
                if lid==4
                    text(Pos_in(lid,1),Pos_in(lid,2),[stars newline opt_axes.Labels{lid} ],'FontSize',12,'HorizontalAlignment','center');
                else
                    text(Pos_in(lid,1),Pos_in(lid,2),[opt_axes.Labels{lid} newline stars],'FontSize',12,'HorizontalAlignment','center');
                end
                
            else
                text(Pos_in(lid,1),Pos_in(lid,2),opt_axes.Labels{lid},'FontSize',12,'HorizontalAlignment','center');
            end
        end
    end
end

uistack(h_iso,'bottom');

hold off;

end

% Function adapt_options adapts the input options struct in a dynamic
% cell-object that can directly pass Name-Value pair arguments through a function
function optList = adapt_options(optStruct)
optList = {};
for optField = fieldnames(optStruct)'
    optList{end+1} = char(optField);                % Name parameter
    optList{end+1} = optStruct.(char(optField));    % Value parameter
end
end