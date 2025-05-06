addpath(genpath("Functions"))

mov = 7;
lw = 2;
legFS = 15;
labFS = 18;
baseFS = 14;
titFS = 20;
cols = colororder;
temp2plot = 'avg';


SOM_ID = '20240318_NCEP_19802010';
load(['SOMs/' SOM_ID '.mat']);
nclasses = DIM_1*DIM_2;
REFYEARS = 1980:2010;

%% Get temperature climatology
tempAll = {'min','avg','max'};
months = monthDays(false);
tBerk=readtable('Data/Processed Reanalysis/SfcAir_InOut/BerkeleyDaily.csv');
dtaRefYear = 1880;
tBerk.WinterYear = tBerk.Year + (tBerk.Month == 12);
Ndays = size(tBerk,1);
climT = zeros(365,3);
tBerkRaw = table2array(tBerk(:,[4 8 5]));
for yr = REFYEARS
    i1 = findDayIndex(yr,1,1,dtaRefYear);
    i2 = findDayIndex(yr,3,1,dtaRefYear);
    i3 = findDayIndex(yr,12,31,dtaRefYear);
    climT(1:59,:) = climT(1:59,:) + tBerkRaw(i1:(i1+58),:);
    climT(60:end,:) = climT(60:end,:) + tBerkRaw(i2:i3,:);
end
climT = climT / length(REFYEARS);
climTmov = runMean(climT,7);
for d = 1:3
    climTmov(d,:) = mean([climT(365-(3-d):365,:); climT(1:d+3,:)],1);
end
for d = 363:365
    climTmov(d,:) = mean([climT(d-3:365,:); climT(1:1-(363-d),:)],1);
end
tAnom = nan(Ndays,3);
for d = 1:Ndays
    if tBerk.Month(d) == 1
        i = tBerk.Day(d);
    elseif tBerk.Month(d) == 2 && tBerk.Day(d) == 29
        i = 59;
    else
        i = sum(months(1:(tBerk.Month(d)-1))) + tBerk.Day(d);
    end
    tAnom(d,:) = tBerkRaw(d,:) - climTmov(i,:);
end
for t = 1:3
    temp = tempAll{t};
    tBerk.(['T' temp 'Anom']) = tAnom(:,t);
end


%% Timeseries Impact
dtaIDAll = {'NCEP','20CRv3','ERA5','ERA20c','JRA55'};
dtaLeg = {'NCEPr1','20CRv3','ERA5','ERA20c','JRA55'};
dtaFields = {'NCEP','CR20','ERA5','ERA20c','JRA55'};
dtaStartYearAll = [1948 1900 1940 1900 1958];
dtaEndYearAll = [2023 2015 2023 2010 2023];

ALL_YEARS = 1881:2022; % berkeley length
NYEARS = length(ALL_YEARS);

Ndta = length(dtaIDAll);
for dta = 1:Ndta
    dtaID = dtaIDAll{dta};
    dtaField = dtaFields{dta};
    TS_IMPACT.(dtaField).min = nan(nclasses,8,NYEARS);
    TS_IMPACT.(dtaField).avg = nan(nclasses,8,NYEARS);
    TS_IMPACT.(dtaField).max = nan(nclasses,8,NYEARS);
    dtaStartYear = dtaStartYearAll(dta);
    dtaEndYear = dtaEndYearAll(dta);
    disp(dtaID)

    NEWSTART = max([dtaStartYearAll(dta)+1 ALL_YEARS(1)]);
    NEWEND = min([dtaEndYearAll(dta) ALL_YEARS(end)]);
    dyears = NEWSTART:NEWEND;
    nyears = length(dyears);
    

    classes = readtable(['Data/Processed Reanalysis/SOM/' SOM_ID '/' dtaID '.csv']);
    classes = classes(classes.WinterYear >= NEWSTART & classes.WinterYear <= NEWEND,:);
    
    tempT = tBerk;
    winter = tempT(((tempT.Month <= 2) | (tempT.Month == 12)) & (tempT.WinterYear >= NEWSTART) & (tempT.WinterYear <= NEWEND),:);
    winterT = table2array(winter(:,11:13));

    YR_REF = classes.WinterYear >= REFYEARS(1) & classes.WinterYear <= REFYEARS(end);
    for Y = 1:nyears
        YR = classes.WinterYear == dyears(Y);
       
        xAll = nan(nclasses,2,3);
        f = nan(nclasses,2);
        for r = 1:2
            if r == 1
                i = YR_REF;
            else
                i = YR;
            end
            for class = 1:nclasses
                classI = i & classes.MapClass == class;
                f(class,r) = sum(classI)/sum(i);
                xAll(class,r,:) = mean(winterT(classI,:),1);
            end
        end
        
        for t = 1:3
            temp = tempAll{t};
            x = xAll(:,:,t);
            thermo = f(:,1) .* diff(x,1,2);
            thermo(isnan(thermo)) = 0;
            dyn = x(:,1) .* diff(f,1,2);
            combined = diff(f,1,2) .* diff(x,1,2);
            combined(isnan(combined)) = 0;
            cum = dyn + thermo + combined;
            ALL = [f x thermo dyn combined cum];
            TS_IMPACT.(dtaField).(temp)(:,:,dyears(Y)-ALL_YEARS(1)+1) = ALL;
        end
    end
end


%% Tmax plot
WH_YEARS = [1958 1988];
TEST_YEARS = [1958 1970 ; 1977 1988];
WH_i = [find(ALL_YEARS == WH_YEARS(1)):find(ALL_YEARS == WH_YEARS(2))];
for t = 1:3
    temp = tempAll{t};
    tssn = seasonFunc(tBerk.(['T' temp 'Anom']),@mean,1881,2022,1880);
    yearOverlap = ALL_YEARS>=REFYEARS(1) & ALL_YEARS<=REFYEARS(end);
    tssn = tssn - mean(tssn(yearOverlap,:),1);
    t_djf.(temp) = tssn(:,1);
    tWH.(temp) = mean(t_djf.(temp)(WH_i));
end

%% Plot
hf = figure;
set(hf, 'Position', [50 50 700 600]);
ha = tight_subplot(3,1,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = ([0.1 0.7 0.82 0.25]);
ha(2).Position = ([0.1 0.4 0.82 0.25]);
ha(3).Position = ([0.1 0.1 0.82 0.25]);
%% Timeseries Plot
for i = 1:3
    axes(ha(i))
    hold on
    xlim([1900 2020])
    box on
    if i ~= 3
        xticklabels([])
        % h = gca;
        % h.XAxis.Visible = 'off';
    end
    ylim([-1.25 1.25])
    yticks(-1:1)
    yticklabels(-1:1)
    set(gca,'FontSize',baseFS,'FontName','Avenir')
    grid on
    yl = ylim();
    xl = TEST_YEARS;
    p = patch(xl(1,[1 2 2 1 1]),yl([1 1 2 2 1]),0.8*[1 1 1],'EdgeColor','none','FaceAlpha',0.25);
    p = patch(xl(2,[1 2 2 1 1]),yl([1 1 2 2 1]),0.8*[1 1 1],'EdgeColor','none','FaceAlpha',0.25);
end

for d = 1:Ndta
    dtaID = dtaFields{d};
    dynTS = squeeze(sum(squeeze(TS_IMPACT.(dtaID).(temp2plot)(:,6,:)),1));
    thermoTS = squeeze(sum(squeeze(TS_IMPACT.(dtaID).(temp2plot)(:,5,:)),1));
    combinedTS = squeeze(sum(squeeze(TS_IMPACT.(dtaID).(temp2plot)(:,7,:)),1));
    
    axes(ha(1))
    plot(ALL_YEARS,runMean(dynTS,mov),"Color",cols(d,:),'LineWidth',lw);
    plot(ALL_YEARS,dynTS,"Color",cols(d,:),'LineWidth',0.25);
    hold on

    axes(ha(2))
    plot(ALL_YEARS,runMean(thermoTS,mov),"Color",cols(d,:),'LineWidth',lw);
    plot(ALL_YEARS,thermoTS,"Color",cols(d,:),'LineWidth',0.25);
    hold on

    axes(ha(3))
    plot(ALL_YEARS,runMean(combinedTS,mov),"Color",cols(d,:),'LineWidth',lw);
    plot(ALL_YEARS,combinedTS,"Color",cols(d,:),'LineWidth',0.25);
    hold on
end

axes(ha(1))
ylabel("Dynamic (°C)",'FontSize',labFS)

axes(ha(2))
ylabel("Thermo (°C)",'FontSize',labFS)

axes(ha(3))
ylabel("Combined (°C)",'FontSize',labFS)
legDta = ["",""];
for d = 1:Ndta
    legDta = [legDta string(dtaLeg{d}) ""];
end
legend(legDta,'NumColumns',6,'Location','southwest','FontSize',legFS)
xlabel("Year",'FontSize',labFS)
xticks([1880:10:2020])
xticklabels([1880:10:2020])
exportgraphics(gcf,['Figures/Figure7.png'],'Resolution',450)
