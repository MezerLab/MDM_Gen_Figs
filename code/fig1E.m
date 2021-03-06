function fig1E(main_dir)

% this function generates Fig.1E

%% Bring Data

% load mixtures MDM measurements:
load(fullfile(main_dir,'/fig1/phantom_data4.mat'));
% This .mat file contains:
% - Slope is a matrix of 12X3 : the rows are different lipid mixtures and
%   the columns are the MDM measurements of each mixture (dR1/dMTV, dR2/dMTV, dMTsat/dMTV).

% load pure samples data:
load(fullfile(main_dir,'/fig1/phantom_data1.mat'));
% This .mat file contains the MRI data of the pure lipids phantoms:
% - LipidMat & LipidMaterror are 55X3X5 matrices contaning mean and std of
%   3 qMRI parameters (see the variable Par) for 5 lipids samples (see the variable lip) with varying water content.
% - MTVMat and MTVMaterror are 55X3X5 matrices contaning mean and std of
%   MTV for 5 lipids samples (see the variable lip) with varying water content.


% The composition of 12 different lipids mixtures (order=PtdCho,PS,Spg):
 w=[0.5, 0  , 0.5 ;...
    2/3, 0  , 1/3 ;...
    1/3, 0  , 2/3 ;...
    0.5, 0.5, 0   ;...
    2/3, 1/3, 0   ;...
    1/3, 2/3, 0   ; ...
    0  , 0.5, 0.5 ;...
    0  , 2/3, 1/3 ;...
    0  , 1/3, 2/3 ;   ...
    1  , 0  , 0   ; ...
    0  ,1   , 0   ;...
    0  ,0   , 1   ];

[n nc]=size(w);
lab={'PtdCho','PS','Spg'}; % lipids in the mixtures
mri_par={'MTsat','R1','R2'}; % qMRI parameters measured

%% set colors

color_earth=[174  182 185; 122 131  144;...
68  82  93; 174  201  141; 104  136  67;...
51  60  30; 241  219  109;...
192  116  6; 104 50  12; 185  153  112;...
105  73  23; 71  40  22;...
152  143  115; 107  91  77; 53  40  31]/255;
colors=color_earth([2,10,5],:);

%% Fig. 1E: Predicting the lipid composition of 12 mixtures using the MDM method. 

% Compute MDM measurements of pure lipids: 
fitobjectT=[];
slopeMat=[];
for jj=1:size(MTVMat,3)
    for ii=1:size(MTVMat,2)
        [v1,I]=sort(MTVMat(:,ii,jj));
        v2=LipidMat(:,ii,jj);
        v2=v2(I);
        v2=v2(~isnan(v1));
        v1=v1(~isnan(v1));
        fitobjectT{ii,jj} = fit(v1,v2,'poly1'); % linear fit qMRI parameters and MTV
        slopeMat(jj,ii)=fitobjectT{ii,jj}.p1;        
    end
end

% Use only the MDM measurements of pure PtdCho, PS, and Spg (as they are in the mixtures): 
Slope_pure=slopeMat([1 3 2],[ 3 1 2]);

% compute cross-validation prediction based on biophysical model: W=[Slope]*[Pure slopes]^-1
prediction_forW=Slope*pinv(Slope_pure);
for ii=1:nc
    mdl = fitlm(w(:,ii),prediction_forW(:,ii));
    R(ii)=mdl.Rsquared.Adjusted;
    P(ii)=mdl.Coefficients{2,4};
    sl(ii)=mdl.Coefficients{2,1};
    int(ii)=mdl.Coefficients{1,1};
end

% generate figure: 
s=figure('Position', [583 489 1110 450]);
for ii=1:nc
    subplot(1,nc,ii)
    hold on
    for r=1:n
        h(r)=scatter(w(r,ii),prediction_forW(r,ii),50,colors(ii,:),'filled');
    end
    hline = refline(1,0);
    hline.Color = 'k';
    hline.LineStyle='--';
    hold off
    xlabel(['true ' lab{ii}]);
    ylabel('prediction');
    hold off
    grid on
    set(gca,'FontSize',12);
    [stars]=pval2stars(P(ii),'num');
    eq={strcat('R^2=',num2str(R(ii),'%0.2f')),[stars]};
    ypos=max(get(gca, 'ylim'));
    xpos=min(get(gca, 'xlim'));
    text(xpos+0.1*xpos,ypos-0.05*ypos,eq,'FontSize',10)
end
set(gcf,'color','w');

end


