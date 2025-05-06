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
set(hf, 'Position', [50 50 1000 450]);
ha = tight_subplot(1,3,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = [0.08 0.15 0.26 0.78];
ha(2).Position = ([0.39 0.15 0.26 0.78]);
ha(3).Position = ([0.7 0.15 0.26 0.78]);

%% Sum Plot
markers = {'v','diamond','^'};
markerSizes = [13 18 13];
spacing = 0.06;
for i = 1:3
    axes(ha(i))
    plot([0 1],[0 0],'k-','LineWidth',1.5)
    xticks([1/6 0.5 5/6])
    xticklabels({'Dynamic','Thermo','Combined'})
    hold on
    ylim([-20 100])
    xlim([0 1])
end
TEST_YEARS = [1958 1970 ; 1977 1988];
TEST_i1 = [find(ALL_YEARS == TEST_YEARS(1,1)):find(ALL_YEARS == TEST_YEARS(1,2))];
TEST_i2 = [find(ALL_YEARS == TEST_YEARS(2,1)):find(ALL_YEARS == TEST_YEARS(2,2))];
for d = 1:5
    dtaID = dtaFields{d};
    dynTS = squeeze(sum(squeeze(TS_IMPACT.(dtaID).(temp2plot)(:,6,:)),1));
    thermoTS = squeeze(sum(squeeze(TS_IMPACT.(dtaID).(temp2plot)(:,5,:)),1));
    combinedTS = squeeze(sum(squeeze(TS_IMPACT.(dtaID).(temp2plot)(:,7,:)),1));

    dynSumTemp.(dtaID) = mean(dynTS(WH_i),'omitnan');
    thermoSumTemp.(dtaID) = mean(thermoTS(WH_i),'omitnan');
    combinedSumTemp.(dtaID) = mean(combinedTS(WH_i),'omitnan');

    dynSum.(dtaID) = 100*mean(dynTS(WH_i),'omitnan') / mean(t_djf.(temp)(WH_i));
    thermoSum.(dtaID) = 100*mean(thermoTS(WH_i),'omitnan') / mean(t_djf.(temp)(WH_i));
    combinedSum.(dtaID) = 100*mean(combinedTS(WH_i),'omitnan') / mean(t_djf.(temp)(WH_i));

    dynSum1.(dtaID) = 100*mean(dynTS(TEST_i1),'omitnan') / mean(t_djf.(temp)(TEST_i1));
    thermoSum1.(dtaID) = 100*mean(thermoTS(TEST_i1),'omitnan') / mean(t_djf.(temp)(TEST_i1));
    combinedSum1.(dtaID) = 100*mean(combinedTS(TEST_i1),'omitnan') / mean(t_djf.(temp)(TEST_i1));

    dynSum2.(dtaID) = 100*mean(dynTS(TEST_i2)) / mean(t_djf.(temp)(TEST_i2));
    thermoSum2.(dtaID) = 100*mean(thermoTS(TEST_i2)) / mean(t_djf.(temp)(TEST_i2));
    combinedSum2.(dtaID) = 100*mean(combinedTS(TEST_i2)) / mean(t_djf.(temp)(TEST_i2));

    axes(ha(1))
    plot([1/6 3/6 5/6]+spacing*(d-3),[dynSum.(dtaID) thermoSum.(dtaID) combinedSum.(dtaID)],'diamond','MarkerSize',markerSizes(2),'MarkerFaceColor',cols(d,:),'MarkerEdgeColor','k');
    axes(ha(2))
    plot([1/6 3/6 5/6]+spacing*(d-3),[dynSum1.(dtaID) thermoSum1.(dtaID) combinedSum1.(dtaID)],'diamond','MarkerSize',markerSizes(2),'MarkerFaceColor',cols(d,:),'MarkerEdgeColor','k');
    axes(ha(3))
    plot([1/6 3/6 5/6]+spacing*(d-3),[dynSum2.(dtaID) thermoSum2.(dtaID) combinedSum2.(dtaID)],'diamond','MarkerSize',markerSizes(2),'MarkerFaceColor',cols(d,:),'MarkerEdgeColor','k');
end
% for t = 1:3
%     temp = tempAll{t};
%     for d = 1:Ndta
%         dtaID = dtaFields{d};
%         plot([1/6 3/6 5/6]+spacing*(d-2.5),[dynSum.(dtaID).(temp) thermoSum.(dtaID).(temp) combinedSum.(dtaID).(temp)],markers{t},'MarkerSize',markerSizes(t),'MarkerFaceColor','none','MarkerEdgeColor','k');
%     end
% end
tits = ["a: WH era (1958-1988)","b: Phase 1 (1958-1970)","c: Phase 2 (1977-1988)"];
for i = 1:3
    axes(ha(i))
    plot([1/3 1/3],ylim(),'LineStyle','-','Color',[0.7500 0.7500 0.7500])
    plot([2/3 2/3],ylim(),'LineStyle','-','Color',[0.7500 0.7500 0.7500])

    %legend(["",string(dtaIDSum)],'Location','northwest','fontsize',legFS)
    set(gca,'FontName','Avenir')
    ax=gca;
    ax.YAxis.FontSize = baseFS;
    ax.XAxis.FontSize = labFS;
    title(tits(i),'FontSize',titFS)
end
axes(ha(1))
ylabel(['WH T_{' temp2plot '} Impact (%)'],'FontSize',labFS)
legend(["",string(dtaLeg)],'Location','northeast','fontsize',legFS)

exportgraphics(gcf,['Figures/Figure6.png'],'Resolution',450)
