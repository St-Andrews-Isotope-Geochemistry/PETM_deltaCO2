cs = load("petmCSMain.mat");

step = .05;
xt = [4:step:9]';
C1 = ksdensity(cs.ClimateSens.values,xt).*step;
%% Figure just Gutjahr
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultTextFontWeight','bold') 

fs = 7;
lw=1;

%colors
c5=[0.831372549019608	0.317647058823529	0.0745098039215686]; %orange
c6=[0.709803921568628	0.698039215686275	0.698039215686275]; %lighter gray

f1 = figure(1); clf;
set(f1,'units','centimeters','pos',[2 10 8.7 5.5],'color','w');
f1.PaperUnits = 'centimeters';
f1.PaperSize = [8.7 5.5];

axes('Position',[.1 .17 .87 .8]);
hold on;
patch([2 2 5 5],[0 0.05 0.05 0],c6,'edgecolor','none');
area(xt,C1,'facecolor',c5,'edgecolor',c5,'FaceAlpha',.5,'linewidth',lw);
set(gca,'fontsize',fs,'box','on','ylim',[0 0.05],'xlim',[4 9]);
text(4.05,.048,'IPCC','fontsize',fs,'horizontalAlignment', 'left');
text(4.05,.044,'AR6 Range','fontsize',fs,'horizontalAlignment', 'left');

xlabel('PETM ECS (˚C per doubling of CO_2)','fontsize',fs);
ylabel('Probability Density','fontsize',fs);

text(7.1,.044,strcat("median ECS = ",sprintf('%0.1f',median(cs.ClimateSens.values)),"˚C"),'FontSize',fs);
text(7.1,.04,strcat("95% CI = ",sprintf('%0.1f',prctile(cs.ClimateSens.values,2.5)),...
    "-",sprintf('%0.1f',prctile(cs.ClimateSens.values,97.5)),"˚C"),'FontSize',fs);

print('ECS_Main.pdf','-painters','-dpdf');