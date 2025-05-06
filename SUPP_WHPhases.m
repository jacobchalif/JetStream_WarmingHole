t=readtable('Data/Processed Reanalysis/SfcAir_InOut/BerkeleyDaily.csv');
Inside_ssn = seasonFunc(t.(['T' 'avg' 'In']),@mean,1881,2022,1880);
years = 1881:2022;
yearRef = years>=1900 & years<=1957;
Inside_anom = Inside_ssn - mean(Inside_ssn(yearRef,:),1);

t=Inside_anom(:,1);
YEARS = 1881:2022;

per = [1901 1957];
avg = mean(t(YEARS >= per(1) & YEARS <= per(2)));
PHASES = [1958 1970 ; 1977 1988];

close all
figure('Position',[39 512 550 350])
plot(YEARS,t,'-*','Color','k')
hold on
plot(YEARS,movmean(t,7),'Color','k','LineWidth',2)
yline(avg,'Color','red','LineWidth',1.1)
xlim([1950 1995])
set(gca,'FontSize',13,'FontName','Avenir')
xlabel("Year","FontSize",18)
ylabel("WH T_{avg} Anomaly (°C)","FontSize",18)
yl = ylim();
xl = PHASES;
p = patch([1958 1988 1988 1958 1958],yl([1 1 2 2 1]),0.8*[1 1 1],'EdgeColor','none','FaceAlpha',0.1,'FaceColor','#0054a8');
p = patch(xl(1,[1 2 2 1 1]),yl([1 1 2 2 1]),0.8*[1 1 1],'EdgeColor','none','FaceAlpha',0.1,'FaceColor','#0054a8');
p = patch(xl(2,[1 2 2 1 1]),yl([1 1 2 2 1]),0.8*[1 1 1],'EdgeColor','none','FaceAlpha',0.1,'FaceColor','#0054a8');
for i = 1:2
    text(mean(PHASES(i,:)),yl(2)-0.5,strcat("Phase ",num2str(i)),'HorizontalAlignment','center','FontSize',18,'FontName','Avenir','FontWeight','bold')
    text(mean(PHASES(i,:)),yl(2)-1,strcat(num2str(PHASES(i,1))," - ",num2str(PHASES(i,2))),'HorizontalAlignment','center','FontSize',13,'FontName','Avenir','FontWeight','bold')
end
text(1964,avg+0.175,"1927-1957 mean",'HorizontalAlignment','center','FontSize',11,'FontName','Avenir','Color','red')
exportgraphics(gcf,'Figures/Supplement/SI5.png','Resolution',450)
