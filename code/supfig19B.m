function supfig19B(diff,diff_M)

% this function generates sup. figure 19B:  The spatial correlation between aging-related changes in volume, iron and water to all MDM dimensions. 

% set color:
c2=[149 143 83;149 170 146;  226 211 112; 222 145 103]./255;
c=[156 163 169; 161 162 120 ; 194 172 124; 152 136 146]./255;

mdlL=[];
pvals=[];

% compute correlations between aging markers:
for ii=2:4
    for jj=2:4
        x=diff(1,:,ii);
        y=diff_M(:,jj);
        mdl = fitlm(x,y);
        % find CI:
        cifig=figure;
        h=plot(mdl);
        CI=[];
        CI(:,1)=h(4).YData;
        CI(:,2)=h(3).YData;
        Xd=h(3).XData;
        close(cifig)
        
        % compute prediction:
        linvec=[min(x),max(x)];
        linePre=linvec.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1}; % fitted line
        linePre2=Xd.*mdl.Coefficients{2,1}+mdl.Coefficients{1,1}; % prediction for each data point
        
        % save fit info:
        boundL{ii,jj}=[linePre2'-CI(:,2),CI(:,1)-linePre2'];
        XdL{ii,jj}=Xd;
        linePre2L{ii,jj}=linePre2;
        linePreL{ii,jj}=linePre;
        linvecL{ii,jj}=linvec;
        mdlL{ii,jj}=mdl;
        pvals(ii,jj)=mdlL{ii,jj}.Coefficients{2,4};
    end
end

% FDR correction for multiple comparisions:
pALLtest=pvals;
pvals_corr=pvals;
pALLtest=pALLtest([2:4],[2:4]);
[FDRval] = mafdr(pALLtest(:),'BHFDR', true);
FDRval=reshape(FDRval,size(pALLtest));
pvals_corr([2:4],[2:4])=FDRval;

% plot figure:
sub=1;
lables={'Chemophysical Effect size','Water Effect size','Volume Effect size','Iron Effect size'};
labelX={'dMTsat/dMTV','dR1/dMTV','dMD/dMTV','dR2/dMTV'};

figure
for ii=2:4
    for jj=2:4
        subplot(3,3,sub)
        hold on
        x=diff(1,:,ii);
        y=diff_M(:,jj);
        [hl,hp]=boundedline(XdL{ii,jj},linePre2L{ii,jj},boundL{ii,jj},'alpha','transparency',0.2,'nan','fill','cmap',[0.5 0.5 0.5],'orientation','vert') ;
        delete(hl);
        scatter(x,y,20,[0.5 0.5 0.5],'filled');
        plot(linvecL{ii,jj},linePreL{ii,jj},'--', 'LineWidth',2,'Color',[0.5 0.5 0.5])
        [stars]=pval2stars(pvals(ii,jj),'num');
        eq=['R^2=',num2str(mdlL{ii,jj}.Rsquared.Adjusted,'%0.2g'),', ', stars];
        xlabel(labelX{ii});
        ylabel(lables{jj});
        title(eq)
        grid on
        sub=sub+1;
        hold off
        ax = gca;
        ax.XColor=c2(ii,:);
        ax.YColor=c(jj,:);
        ax.FontWeight='bold';
    end
end
set(gcf, 'Position',[1 1 1053 800]);

end