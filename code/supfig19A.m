
function supfig19A(diff_M,diff)

% this function generates sup. figure 19A: Distinct spatial patterns of different aging markers throughout the brain

% set order of ROIs:
ord=1:14;
ord([5:7])=[1,3,4];
ord([1:4])=[11:14];
ord([8:14])=[5:6,2,9:10,7:8];
[M,I] = max(abs(diff_M(ord,:)),[],2);

% set color:
c=[156 163 169; 161 162 120 ; 194 172 124; 152 136 146]./255;
c2=[149 143 83;149 170 146;  226 211 112; 222 145 103]./255;

% plot figure:
figure
hold on
for ii=1:14
    rectangle('Position',[(ii)-0.5 -2.5 1 8],'FaceColor',[c(I(ii),:) 0.5],'EdgeColor','none')
end
rectangle('Position',[1-0.75 -0.2 14.76 0.4],'FaceColor',[1 1 1 0.7],'LineWidth',0.005,'EdgeColor','k','LineStyle','--')

plot([0 14.5],[0 0],'k')

for ii=1:4
    s(ii)=plot((diff_M(ord,ii)),'-*','lineWidth',4,'Color',c(ii,:));
end

for ii=2:4
    s(3+ii)=plot(diff(1,ord,ii),'-*','lineWidth',4,'Color',c2(ii,:));
end
hold off
set(gca,'FontSize',15)
legend(s,{'dMTsat/dMTV','Water','Volume','Iron','dR1/dMTV','dMD/dMTV','dR2/dMTV'},'Color',[1 1 1 0.5],'EdgeColor','none','FontSize',15)
set(gca, 'xtick',1:14, 'xticklabel', '','FontSize',15);
ylabel('Effect size','fontSize',18);
xlim([0.5 14.5]);
ax = gca;
ax.YGrid = 'on';
ylim([-2 3.5]);
set(gcf, 'Position',[1 1 1053 618]);




end