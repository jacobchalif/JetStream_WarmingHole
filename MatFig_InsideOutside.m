addpath(genpath("Functions/"))

temp = 'avg'; % min, avg, max
mov = 7;
lw = 1.25;
YL = [-3.5 3.5];
titFS = 20;
legFS = 16;
labFS = 18;
baseFS = 14;
CSTEP = 0.25;
CMIN = -2;
CMAX = 2;
CS = [CMIN:CSTEP:(-1*CSTEP) CSTEP:CSTEP:CMAX];

%% Plotting
dtasets = {'NCEP','20CRv3','ERA5','ERA20c','JRA55'};
dtasetLabels = {'NCEPr1','20CRv3','ERA5','ERA20c','JRA55'};
starts = [1948 1900 1940 1900 1958];
ends = [2023 2015 2023 2010 2023];
cols = colororder();
cols(end+1,:) = [0 0 0];
overlap = [max(starts) min(ends)];
overlap = [1980 2010];
legDta = {};
for i = 1:length(dtasets)
    legDta{end+1} = dtasets{i};
    legDta{end+1} = dtasets{i};
end
ssn = 1; % 1 = DJF

hf = figure;
set(hf, 'Position', [50 50 1200 900]);
ha = tight_subplot(2,2,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = ([0.046 0.52 0.55 0.44]);
ha(2).Position = ([0.65 0.6 0.285 0.3]);
ha(3).Position = ([0.2 0.1 0.7 0.3]);
ha(4).Position = ([0 0 0 0]);

%% Panel A
lat = readmatrix('Data/Processed Reanalysis/Gridded/Berkeley Daily/lat.csv');
lon = readmatrix('Data/Processed Reanalysis/Gridded/Berkeley Daily/lon.csv');
tAnom = readmatrix(['Data/Processed Reanalysis/Gridded/Berkeley Daily/T' temp '_1958_1988_base_1900_1957.csv']);

axes(ha(1))
axis off;
ax=usamap('conus');
ax.Position= ha(1).Position;
axes(ax)
states = shaperead('usastatelo', 'UseGeoCoords', true);
%p=pcolorm(lat,lon,tAnom);
%p.EdgeColor='none';
contourfm(lat,lon,tAnom,CS);
geoshow(states,'FaceColor',[1,1,1],'facealpha',0); %Adds state boundarie
colormap(ax,redblue)
caxis([-2 2])


% draws rectangle around WH
patchm([30 40 40 30 30],[-95 -95 -82 -82 -95],'LineWidth',3,'LineStyle','--','FaceColor','none')

% % X axis information
% xlim([-125 -65]);
% xticks(-120:10:-70)
% xticklabels(num2str((-120:10:-70)'))
% 
% % Y axis informatio
% ylim([25 50]);   
% yticks(25:5:50)
% yticklabels(num2str((25:5:50)'))
A = colorbar('southoutside'); % put colorbar at the bottom
set(A, 'Position', [0.14 .52 .365 .015]); %[left,bottom,width,height] [0.33 .08 .365 .015]
A.FontSize = baseFS+1;
A.FontName='Avenir';
A.Label.String = ['DJF T_{' temp '}, 1958-1988 minus 1927-1957 (°C)']; % name colorbar title
A.Label.FontSize = labFS; 
A.Label.FontName = 'Avenir'; 

a=annotation('textbox',[0.13 .94 .365 .015],'String',['a: Winter T_{' temp '} difference'],'FontName','Avenir');
a.FontSize = titFS;
a.EdgeColor = 'none';
a.HorizontalAlignment = 'center';
a.FontWeight = 'bold';

%% Panel B:
t=readtable('Data/Processed Reanalysis/SfcAir_InOut/BerkeleyDaily.csv');
Inside_ssn = seasonFunc(t.(['T' temp 'In']),@mean,1881,2022,1880);
Outside_ssn = seasonFunc(t.(['T' temp 'Out']),@mean,1881,2022,1880);
years = 1881:2022;
yearRef = years>=1927 & years<=1957;
Inside_anom = Inside_ssn - mean(Inside_ssn(yearRef,:),1);
Outside_anom = Outside_ssn - mean(Outside_ssn(yearRef,:),1);
    
axes(ha(2))
set(gca,'FontSize',baseFS,'FontName','Avenir')
grid on
box on

y21=Inside_anom(:,ssn); % tmax inside WH
y22=Outside_anom(:,ssn); % tmax outside WH

yy21 = runMean(y21,7); % smooth inside WH tmax
yy22 = runMean(y22,7); % smooth outside WH tmax

% plots unsmoothed data
patchline(years,y21,'edgecolor','blue','linewidth',1,'linestyle','--','edgealpha',0.1);
hold on
patchline(years,y22,'edgecolor','red','linewidth',1,'linestyle','--','edgealpha',0.1);

plot(years,yy21,'r-','Color','blue','LineWidth',1.5);
plot(years,yy22,'r-','Color','red','LineWidth',1.5);

% X axis information
xlim([1900 2023]);
xticks(1900:20:2020)
xticklabels(num2str((1900:20:2020)'))
xlabel('Year','fontsize',labFS);

% Y axis informatio
ylim(YL)
yticks(-3:3)
yticklabels(num2str((-3:3)'))
ylabel(['T_{' temp '} Anomaly (°C)'],'fontsize',labFS);


title(['b: Winter T_{' temp '}'], 'FontSize',titFS);  
lgd = legend('','',['T_{' temp '} inside WH'],['T_{' temp '} outside WH']);
lgd.FontSize = legFS;
lgd.Location='northwest';
legend boxoff


%% Panels C and D:
axes(ha(3))
for d = 1:length(dtasets)
    disp(dtasets{d})
    t=readtable(['Data/Processed Reanalysis/SfcAir_InOut/' dtasets{d} '.csv']);
    Inside_ssn = seasonFunc(t.(['T' temp 'In']),@mean,starts(d),ends(d),starts(d));
    Outside_ssn = seasonFunc(t.(['T' temp 'Out']),@mean,starts(d),ends(d),starts(d));
    years = starts(d):ends(d);
    yearOverlap = years>=overlap(1) & years<=overlap(2);
    Inside_anom = Inside_ssn - mean(Inside_ssn(yearOverlap,:),1);
    Outside_anom = Outside_ssn - mean(Outside_ssn(yearOverlap,:),1);
    Diff_anom = Inside_anom - Outside_anom;

    axes(ha(3))
    plot(years,runMean(Inside_anom(:,ssn),mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
end

if string(temp) == "avg"
    b = readtable(['Data/Processed Reanalysis/SfcAir_InOut/GISTEMP.csv']);
    years = 1881:2023;
    tIn = nan(length(years),1);
    tOut = nan(length(years),1);
    for yr = years
        i = (b.Year == yr & b.Month <= 2) | (b.Year == (yr-1) & b.Month == 12);
        tIn(yr-years(1)+1) = mean(b.Tin(i));
        tOut(yr-years(1)+1) = mean(b.Tout(i));
    end
    yearOverlap = years>=overlap(1) & years<=overlap(2);
    tIn_anom = tIn - mean(tIn(yearOverlap));
    tOut_anom = tOut - mean(tOut(yearOverlap));
    axes(ha(3))
    plot(years,runMean(tIn_anom,7),'k','LineWidth',2','Color','k','LineStyle','--')
 
    dtasetLabels{end+1} = 'GISTEMP4';
    
    
    b = readtable(['Data/Processed Reanalysis/SfcAir_InOut/CRUTEM5.csv']);
    years = 1851:2023;
    tIn = nan(length(years),1);
    tOut = nan(length(years),1);
    for yr = years
        i = (b.Year == yr & b.Month <= 2) | (b.Year == (yr-1) & b.Month == 12);
        tIn(yr-years(1)+1) = mean(b.Tin(i));
        tOut(yr-years(1)+1) = mean(b.Tout(i));
    end
    yearOverlap = years>=overlap(1) & years<=overlap(2);
    tIn_anom = tIn - mean(tIn(yearOverlap));
    tOut_anom = tOut - mean(tOut(yearOverlap));
    axes(ha(3))
    plot(years,runMean(tIn_anom,7),'k','LineWidth',2','Color','k','LineStyle',':')
    
    dtasetLabels{end+1} = 'CRUTEM5';
end

b = readtable(['Data/Processed Reanalysis/SfcAir_InOut/USHCN.FLs.52j.csv']);
years = 1900:2023;
tIn = nan(length(years),1);
tOut = nan(length(years),1);
for yr = years
    i = (b.Year == yr & b.Month <= 2) | (b.Year == (yr-1) & b.Month == 12);
    tIn(yr-years(1)+1) = mean(b.(['T' temp 'In'])(i));
    tOut(yr-years(1)+1) = mean(b.(['T' temp 'Out'])(i));
end
yearOverlap = years>=overlap(1) & years<=overlap(2);
tIn_anom = tIn - mean(tIn(yearOverlap));
tOut_anom = tOut - mean(tOut(yearOverlap));

axes(ha(3))
plot(years,runMean(tIn_anom,7),'LineWidth',lw','Color','#949494','LineStyle','-')


dtasetLabels{end+1} = 'USHCN';
   
b = readtable(['Data/Processed Reanalysis/SfcAir_InOut/BerkeleyDaily.csv']);
years = 1851:2023;
tIn = nan(length(years),1);
tOut = nan(length(years),1);
for yr = years
    i = (b.Year == yr & b.Month <= 2) | (b.Year == (yr-1) & b.Month == 12);
    tIn(yr-years(1)+1) = mean(b.(['T' temp 'In'])(i));
    tOut(yr-years(1)+1) = mean(b.(['T' temp 'Out'])(i));
end
yearOverlap = years>=overlap(1) & years<=overlap(2);
tIn_anom = tIn - mean(tIn(yearOverlap));
tOut_anom = tOut - mean(tOut(yearOverlap));

axes(ha(3))
plot(years,runMean(tIn_anom,7),'k','LineWidth',2','Color','r','LineStyle','--')
    
dtasetLabels{end+1} = 'Berkeley Earth (daily)';


i2 = 3;
for i = 3:i2
    axes(ha(i))
    set(gca,'FontSize',baseFS,'FontName','Avenir')
    xlim([1900 2023])
    ylim(YL)
    lgd = legend(dtasetLabels,'Location','eastoutside');
    lgd.FontSize = legFS;    
    grid on
end

axes(ha(3))
ylabel(['T_{' temp '} Anomaly (°C)'],'fontsize',labFS);
text(1905,2.8,['c: Winter T_{' temp '} inside WH'],'fontsize',titFS,'fontweight','bold','FontName','Avenir')
xlabel("Year",'fontsize',labFS);

%% LINES
annotation('line',[0.43 0.65],[0.678 0.6],'LineWidth',1.5,'LineStyle',':')
annotation('line',[0.414 0.65],[0.793 0.9],'LineWidth',1.5,'LineStyle',':')

exportgraphics(gcf,['Figures/Figure2.png'],'Resolution',450)


figure('Position',[100 100 700 70])
c=colormap('redblue');
%CS = [-1.2:0.2:-0.2 0.2:0.2:1.2];
NC = length(CS)-1;

cols = nan(NC,3);
for i = 1:NC
   x1 = CS(i);
   x2 = CS(i+1);
   cols(i,:) = c(1+floor((i-1)*(length(c)-1)/(NC-1)),:);
   patch([x1 x1 x2 x2 x1],[0 1 1 0 0],cols(i,:))
end
xlim([CS(2) CS(end-1)])
xticks(CS)
yticks([])
set(gca,'FontSize',baseFS,'FontName','Avenir')

exportgraphics(gcf,['Figures/Colorbars/Figure2.png'],'Resolution',450)
