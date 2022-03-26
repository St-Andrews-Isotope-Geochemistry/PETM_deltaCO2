csP = load("petmCSMain.mat");
cs = load("petmCSAlk.mat");
csO = load("petmCSOmega.mat");

%ECS
step = .05;
xt = [3:step:11]';
O1 = ksdensity(csO.ClimateSens.values,xt).*step;
P1 = ksdensity(csP.ClimateSens.values,xt).*step;
C1 = ksdensity(cs.ClimateSens.values,xt).*step;

%doublings
step = .005;
xt2 = [0.45:step:2]';
O2 = ksdensity(csO.combined_doublings.sampler.samples,xt2).*step;
P2 = ksdensity(csP.combined_doublings.sampler.samples,xt2).*step;
C2 = ksdensity(cs.combined_doublings.sampler.samples,xt2).*step;

%deltaCO2
step2 = 20;
xt3 = [100:step2:2300]';
O3 = ksdensity(csO.combined_delCO2.sampler.samples,xt3).*step2;
P3 = ksdensity(csP.combined_delCO2.sampler.samples,xt3).*step2;
C3 = ksdensity(cs.combined_delCO2.sampler.samples,xt3).*step2;

%% Figure
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultTextFontWeight','bold') 

fs = 7;
lw=1;

%colors
c1=hex2rgb('#95190C'); %red
c2=hex2rgb('#E28413'); %yellow
c5=hex2rgb('#D45113'); %orange
c6=hex2rgb('#B5B2B2'); %lighter gray

f1 = figure(1); clf;
set(f1,'units','centimeters','pos',[2 10 17.8 5.5],'color','w');
f1.PaperUnits = 'centimeters';
f1.PaperSize = [17.8 5.5];

wd = 0.85./3;
ht = 0.83;
s1 = 0.05;
v1 = 0.14;
hoff = 0.04;

%ECS

axes('Position',[s1 v1 wd ht]);
hold on;
patch([2 2 5 5],[0 0.05 0.05 0],c6,'edgecolor','none');
a1=area(xt,O1,'facecolor',c2,'edgecolor',c2,'FaceAlpha',.5,'linewidth',lw);
a3=area(xt,C1,'facecolor',c1,'edgecolor',c1,'FaceAlpha',.5,'linewidth',lw);
a2=area(xt,P1,'facecolor',c5,'edgecolor',c5,'FaceAlpha',.5,'linewidth',lw);
set(gca,'fontsize',fs,'box','on','ylim',[0 0.05],'xlim',[3.5 10.5]);
L1=legend([a3 a2 a1],'Alk','Main','\Omega','box','off');
L1.ItemTokenSize = [8 18];
L1.FontSize = fs;
L1.Position = [.25 .8 .1 .1];
text(3.7,.041,'IPCC','fontsize',fs,'horizontalAlignment', 'left');
text(3.7,.038,'AR6','fontsize',fs,'horizontalAlignment', 'left');
text(3.7,.048,'a.','fontsize',fs);

xlabel('PETM ECS (^{\circ}C per doubling of CO_2)','fontsize',fs);
yL = ylabel('Probability Density','fontsize',fs);

%Doublings
axes('Position',[s1+wd+hoff v1 wd ht]);
hold on;
a1=area(xt2,O2,'facecolor',c2,'edgecolor',c2,'FaceAlpha',.5,'linewidth',lw);
a2=area(xt2,P2,'facecolor',c5,'edgecolor',c5,'FaceAlpha',.5,'linewidth',lw);
a3=area(xt2,C2,'facecolor',c1,'edgecolor',c1,'FaceAlpha',.5,'linewidth',lw);
set(gca,'fontsize',fs,'box','on','xlim',[0.5 1.3],'xtick',[.5:.25:1.5]);
L1=legend([a3 a2 a1],'Alk','Main','\Omega','box','off');
L1.ItemTokenSize = [8 18];
L1.FontSize = fs;
L1.Position = [.25+wd+hoff .8 .1 .1];
xlabel('CO_2 doublings','fontsize',fs);
yL = ylabel('Probability','fontsize',fs);
yL.Position = [.49 3 -1];
text(.52,.038,'b.','fontsize',fs);

%prePETM
axes('Position',[s1+(wd+hoff).*2 v1 wd ht]);
hold on;
a1=area(xt3,O3,'facecolor',c2,'edgecolor',c2,'FaceAlpha',.5,'linewidth',lw);
a2=area(xt3,P3,'facecolor',c5,'edgecolor',c5,'FaceAlpha',.5,'linewidth',lw);
a3=area(xt3,C3,'facecolor',c1,'edgecolor',c1,'FaceAlpha',.5,'linewidth',lw);
set(gca,'fontsize',fs,'box','on','xlim',[200 2000]);
L1=legend([a3 a2 a1],'Alk','Main','\Omega','box','off');
L1.ItemTokenSize = [8 18];
L1.FontSize = fs;
L1.Position = [.25+(wd+hoff).*2 .8 .1 .1];
xlabel('\Delta CO_{2}','fontsize',fs);
text(250,.067,'c.','fontsize',fs);

print('doublings.pdf','-painters','-dpdf');