addpath(genpath("Functions/"))

option = 2;
    

mov = 7;
lw = 1.25;
YL = [-3.5 3.5];
titFS = 18;
legFS = 14;
labFS = 16;
baseFS = 12;
        

SOM_ID = '20240318_NCEP_19802010';
load(['SOMs/' SOM_ID '.mat']);
nclasses = DIM_1*DIM_2;

dtasets = {'NCEP','20CRv3','ERA5','ERA20c','JRA55'};
dtasetLabels = {'NCEPr1','20CRv3','ERA5','ERA20c','JRA55'};
legDta = dtasetLabels;
starts = [1948 1900 1940 1900 1958];
ends = [2023 2015 2023 2010 2023];
cols = colororder();
cols(end+1,:) = [0 0 0];
overlap = [max(starts) min(ends)];
overlap = [1980 2010];
ssn = 1; % 1 = DJF

t=readtable('Data/Processed Reanalysis/SfcAir_InOut/BerkeleyDaily.csv');
Inside_ssn = seasonFunc(t.TavgIn,@mean,1881,2022,1880);
Tyears = 1881:2022;
yearOverlap = Tyears>=1980 & Tyears<=2010;
T_anom = Inside_ssn - mean(Inside_ssn(yearOverlap,:),1);

%% Selecting wavy classes
LWA = readtable("Data/Processed Reanalysis/LWA/NCEP.csv");
MCI = readtable("Data/Processed Reanalysis/MCI300/NCEP.csv");
classes = readmatrix(['Data/Processed Reanalysis/SOM/' SOM_ID '/master.csv']);
N = findDayIndex(2023,12,31,1948);
dates = nan(N,3);
for year = 1948:2023
    l = isLeap(year);
    m = monthDays(l);
    i1 = findDayIndex(year,1,1,1948);
    i2 = findDayIndex(year,12,31,1948);
    years(i1:i2) = year;
    dates(i1:i2,1) = year;
    for i = 1:12
        if i == 1
            dates(i1:(i1+m(i)-1),2) = i;
            dates(i1:(i1+m(i)-1),3) = 1:m(i);
        else
            dates((i1+sum(m(1:(i-1)))):(i1+sum(m(1:(i-1)))+m(i)-1),2) = i;
            dates((i1+sum(m(1:(i-1)))):(i1+sum(m(1:(i-1)))+m(i)-1),3) = 1:m(i);
        end
    end
end
winterI = ((dates(:,1) >= 1980) & (dates(:,1) <= 2010)) & ((dates(:,2) == 12) | (dates(:,2) <= 2));
dates = dates(winterI,:);
LWA = LWA(winterI,:);
MCI = MCI(winterI,:);
threshold = 0.8;
classWaves = false(nclasses,1);
for class = 1:nclasses
    classI = classes(:,4) == class;
    n = sum(classI);
    LWA_class = LWA.LWA_Anticylonic(classI) - LWA.LWA_Cylonic(classI);
    MCI_class = MCI.MCI(classI);
    if ((sum(LWA_class<0)/n) >= threshold) && ((sum(MCI_class<0)/n) >= threshold)
        classWaves(class) = true;
    end
end
WAVES = find(classWaves);


%% just timeseries
hf = figure;
set(hf, 'Position', [50 50 500 1000]);
ha = tight_subplot(4,1,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = ([0.18 0.73 0.66 0.18]);
ha(2).Position = ([0.18 0.52 0.66 0.18]);
ha(3).Position = ([0.18 0.31 0.66 0.18]);
ha(4).Position = ([0.18 0.1 0.66 0.18]);


corrs = nan(length(dtasets),3);
ps = nan(length(dtasets),3);
for d = 1:length(dtasets)
    disp(dtasets{d})
    years = starts(d):ends(d);
    yearOverlap = years>=overlap(1) & years<=overlap(2);

    %% calculating waviness frequency
    C=readtable(['Data/Processed Reanalysis/SOM/' SOM_ID '/' dtasets{d} '.csv']);
    classes = C.MapClass;
    classYrs = nan(length(years),nclasses);
    for class = 1:nclasses
        classI = classes == class;
        for yr = 1:length(years)
            yrI = C.WinterYear == years(yr);
            classYrs(yr,class) = sum(classI & yrI) / sum(yrI);
        end
    end
    waves = sum(classYrs(:,WAVES),2);
    [corrs(d,1), ps(d,1)] = corr(waves(years <= 2022),T_anom(Tyears >= years(1) & Tyears <= years(end))');

    % panel a: troughing
    axes(ha(1))
    plot(years,runMean(waves,mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir','YDir','reverse')
    
    % panel b: LWA
    lwa=readtable(['Data/Processed Reanalysis/LWA/' dtasets{d} '.csv']);
    lwa_ssn = seasonFunc(lwa.LWA_Anticylonic-lwa.LWA_Cylonic,@mean,starts(d),ends(d),starts(d));
    lwa_anom = lwa_ssn - mean(lwa_ssn(yearOverlap,:),1);
    [corrs(d,2), ps(d,2)] = corr(lwa_anom(years <= 2022,1),T_anom(Tyears >= years(1) & Tyears <= years(end))');

    axes(ha(2))
    plot(years,runMean(lwa_anom(:,ssn)/1e6,mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir')

    % panel c: MCI
    mci=readtable(['Data/Processed Reanalysis/MCI300/' dtasets{d} '.csv']);
    mci_ssn = seasonFunc(mci.MCI,@mean,starts(d),ends(d),starts(d));
    mci_anom = mci_ssn - mean(mci_ssn(yearOverlap,:),1);
    [corrs(d,3), ps(d,3)] = corr(mci_anom(years <= 2022,1),T_anom(Tyears >= years(1) & Tyears <= years(end))');

    axes(ha(3))
    plot(years,runMean(mci_anom(:,ssn),mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir')

    % panel 3: zonal winds
    uwnd=readtable(['Data/Processed Reanalysis/U300/' dtasets{d} '.csv']);
    uwnd_ssn = seasonFunc(uwnd.UWnd,@mean,starts(d),ends(d),starts(d));
    uwnd_anom = uwnd_ssn - mean(uwnd_ssn(yearOverlap,:),1);
    [corrs(d,3), ps(d,3)] = corr(uwnd_anom(years <= 2022,1),T_anom(Tyears >= years(1) & Tyears <= years(end))');

    axes(ha(4))
    plot(years,runMean(uwnd_anom(:,ssn),mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir')
    xlabel("Year",'FontSize',labFS)
end
for i = 1:4
    axes(ha(i))
    xticks(1900:20:2020)
    xlim([1900 2020])
    grid on
    if i < 4
        xticklabels([])
    end
end
%% temperature
legDta{end+1} = 'T_{avg}';

for i = 1:4
    axes(ha(i))
    yyaxis right
    plot(Tyears,movmean(T_anom(:,1),mov),'LineWidth',lw,'Color','k')
    set(gca,'YColor','k')
    ylabel("DJF T_{avg} (°C)",'FontSize',labFS)
    ylim([-2 2])
end
ty = 1.65;
tx = 1905;

axes(ha(1))
yyaxis left
ylabel("DJF STI",'FontSize',labFS)
ylim([0.1 0.5])
yticks(0.1:0.1:0.5)
legend(legDta,'NumColumns',4,'Location','northoutside','FontSize',legFS)
yyaxis right
text(tx,ty,"a: SOM Troughing Index",'FontSize',titFS,'FontName','Avenir')

axes(ha(2))
yyaxis left
ylabel("DJF LWA anom. (10^{6} m^{2})",'FontSize',labFS)
ylim([-5 5])
yticks(-5:2.5:5)
ha(1).Position = ([0.15 0.735 0.72 0.19]);
ha(2).Position = ([0.15 0.505 0.72 0.19]);
ha(3).Position = ([0.15 0.275 0.72 0.19]);
ha(4).Position = ([0.15 0.045 0.72 0.19]);

yyaxis right
text(tx,ty,"b: LWA",'FontSize',titFS,'FontName','Avenir')

axes(ha(3))
yyaxis left
ylabel("DJF MCI anom.",'FontSize',labFS)
yyaxis right
text(tx,ty,"c: MCI",'FontSize',titFS,'FontName','Avenir')


axes(ha(4))
yyaxis left
ylabel("DJF U-wind anom. (m/s)",'FontSize',labFS)
yyaxis right
text(tx,ty,"d: U-wind",'FontSize',titFS,'FontName','Avenir')

exportgraphics(gcf,['Figures/Figure5.png'],'Resolution',450)


