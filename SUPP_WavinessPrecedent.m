addpath(genpath("Functions/"))    

mov = 1;
lw = 0.75;
YL = [-3.5 3.5];
titFS = 18;
legFS = 14;
labFS = 16;
baseFS = 12;
        

SOM_ID = '20240318_NCEP_19802010';
load(['SOMs/' SOM_ID '.mat']);
nclasses = DIM_1*DIM_2;

dtasets = {'NCEP','20CRv3','ERA5','ERA20c','JRA55'};%,'MERRA2','NCEP2'};
dtasetLabels = {'NCEPr1','20CRv3','ERA5','ERA20c','JRA55'};%,'MERRA2','NCEP2'};
legDta = dtasetLabels;
starts = [1948 1900 1940 1900 1958];% 1980 1979];
ends = [2023 2015 2023 2010 2023];% 2023 2023];
cols = colororder();
cols(end+1,:) = [0 0 0];
overlap = [max(starts) min(ends)];
overlap = [1980 2010];
ssn = 1; % 1 = DJF

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
set(hf, 'Position', [50 50 500 800]);
ha = tight_subplot(3,1,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = ([0.18 0.62 0.66 0.25]);
ha(2).Position = ([0.18 0.35 0.66 0.25]);
ha(3).Position = ([0.18 0.08 0.66 0.25]);

minimaYrs = [2009 2010];

for d = 1:length(dtasets)
    disp(dtasets{d})
    years = starts(d):ends(d);
    yearOverlap = years>=overlap(1) & years<=overlap(2);
    minI = ismember(years,minimaYrs);
    minLW = 2;

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

    % panel a: troughing
    axes(ha(1))
    plot(years,runMean(waves,mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir','YDir','reverse')
    minima = mean(waves(minI));
    yline(minima,'LineWidth',minLW,'Color',cols(d,:),'LineStyle','--')
    below = waves > minima;
    scatter(years(below),waves(below),40,cols(d,:),'filled','o')
    numBelow = sum(below(years >= 1960 & years < 2000));
    annotation( 'textbox', 'String', strcat(num2str(numBelow)," (",num2str(round(100*numBelow/40))," %)"), ...
            'Color', cols(d,:), 'FontSize', 14, 'FontName','Avenir', ...
            'EdgeColor', 'none', 'Position', [0.85,0.685+d/40,0.15,0] )
    
    % panel b: LWA
    lwa=readtable(['Data/Processed Reanalysis/LWA/' dtasets{d} '.csv']);
    lwa_ssn = seasonFunc(lwa.LWA_Anticylonic-lwa.LWA_Cylonic,@mean,starts(d),ends(d),starts(d));
    lwa_anom = lwa_ssn - mean(lwa_ssn(yearOverlap,:),1);

    axes(ha(2))
    plot(years,runMean(lwa_anom(:,ssn)/1e6,mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir')
    minima = mean(lwa_anom(minI,ssn)/1e6);
    yline(minima,'LineWidth',minLW,'Color',cols(d,:),'LineStyle','--')
    below = lwa_anom(:,ssn)/1e6 < minima;
    scatter(years(below),lwa_anom(below,ssn)/1e6,40,cols(d,:),'filled','o')
    numBelow = sum(below(years >= 1960 & years < 2000));
    annotation( 'textbox', 'String', strcat(num2str(numBelow)," (",num2str(round(100*numBelow/40))," %)"), ...
            'Color', cols(d,:), 'FontSize', 14, 'FontName','Avenir', ...
            'EdgeColor', 'none', 'Position', [0.85,0.41+d/40,0.15,0] )

    % panel c: MCI
    mci=readtable(['Data/Processed Reanalysis/MCI300/' dtasets{d} '.csv']);
    mci_ssn = seasonFunc(mci.MCI,@mean,starts(d),ends(d),starts(d));
    mci_anom = mci_ssn - mean(mci_ssn(yearOverlap,:),1);

    axes(ha(3))
    plot(years,runMean(mci_anom(:,ssn),mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir')
    minima = mean(mci_anom(minI,ssn));
    yline(minima,'LineWidth',minLW,'Color',cols(d,:),'LineStyle','--')
    below = mci_anom(:,ssn) < minima;
    scatter(years(below),mci_anom(below,ssn),40,cols(d,:),'filled','o')
    numBelow = sum(below(years >= 1960 & years < 2000));
    annotation( 'textbox', 'String', strcat(num2str(numBelow)," (",num2str(round(100*numBelow/40))," %)"), ...
            'Color', cols(d,:), 'FontSize', 14, 'FontName','Avenir', ...
            'EdgeColor', 'none', 'Position', [0.85,0.145+d/40,0.15,0] )

    xlabel("Year",'FontSize',labFS)
end
for i = 1:3
    axes(ha(i))
    xticks(1900:10:2020)
    xlim([1960 2020])
    grid on
    if i < 3
        xticklabels([])
    end
end
tx = 1962;

axes(ha(1))
ylabel("STI",'FontSize',labFS)
ylim([0 0.6])
yticks(0:0.1:0.6)
legDta2 = {};
for i = 1:length(legDta)
    legDta2{end+1} = legDta{i};
    legDta2{end+1} = '';
    legDta2{end+1} = '';
end
legend(legDta2,'NumColumns',4,'Location','northoutside','FontSize',legFS)
text(tx,0.05,"a: SOM Troughing Index",'FontSize',titFS,'FontName','Avenir')

axes(ha(2))
ylabel("LWA anomaly (10^{6} m^{2})",'FontSize',labFS)
ylim([-10 7.5])
yticks(-10:5:7.5)
ha(1).Position = ([0.12 0.62 0.72 0.25]);
ha(2).Position = ([0.12 0.35 0.72 0.25]);
ha(3).Position = ([0.12 0.08 0.72 0.25]);
text(tx,6,"b: LWA",'FontSize',titFS,'FontName','Avenir')

axes(ha(3))
ylabel("MCI anomaly",'FontSize',labFS)
text(tx,0.125,"c: MCI",'FontSize',titFS,'FontName','Avenir')

exportgraphics(gcf,'Figures/Supplement/SI7.png','Resolution',450)
