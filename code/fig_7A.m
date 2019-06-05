
function fig_6A(diff_M,std_M)

% this function generates figure 7A: Distinct spatial patterns of different aging markers throughout the brain

% set order of ROIs
ord=1:14;
ord([5:7])=[1,3,4];
ord([1:4])=[11:14];
ord([8:14])=[5:6,2,9:10,7:8];
[M,I] = max(abs(diff_M(ord,:)),[],2);

% color:
c=[156 163 169; 161 162 120 ; 194 172 124; 152 136 146]./255;

% plot figure:
figure(1);
hold on
for ii=1:14
    rectangle('Position',[(ii)-0.5 -2.5 1 8],'FaceColor',[c(I(ii),:) 0.5],'EdgeColor','none')
end
rectangle('Position',[1-0.75 -0.2 14.76 0.4],'FaceColor',[1 1 1 0.7],'LineWidth',0.005,'EdgeColor','k','LineStyle','--')
plot([0 14.5],[0 0],'k')
for ii=1:4
    s(ii)=plot((diff_M(ord,ii)),'-*','lineWidth',4,'Color',c(ii,:));
    [hl,hp]=boundedline(1:14,diff_M(ord,ii),std_M(ord,ii),'alpha','transparency',0.5,'nan','fill','cmap',c(ii,:),'orientation','vert') ;
end
hold off
set(gca,'FontSize',15)
legend(s,{'Chemophysical','Water','Atrophy','Iron'},'Color',[1 1 1 0.5],'EdgeColor','none','FontSize',18)
set(gca, 'xtick',1:14, 'xticklabel', '','FontSize',15);
ylabel('Effect size','fontSize',18);
xlim([0.5 14.5]);
ax = gca;
ax.YGrid = 'on';
ylim([-2.5 2.5]);
set(gcf, 'Position',[1 1 1053 618]);


end