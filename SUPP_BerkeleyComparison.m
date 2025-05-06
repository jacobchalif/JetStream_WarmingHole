addpath(genpath("Functions"))

REF = [1980 2010];

dtasets = {'Berkeley','CRUTEM5','GISTEMP','USHCN.FLs.52j','BerkeleyDaily','Livneh','USHCN.raw','BerkeleyNotQC'};
dtasetNames = {'BerkeleyMonthly','CRUTEM5','GISTEMP4','USHCNcorr','BerkeleyDaily','Livneh','USHCNraw','BerkeleyRaw'};
dtasetLegend = {'Berkeley (m, c)','CRUTEM5 (m, c)','GISTEMP4 (m, c)','USHCN (m, c)','Berkeley (d, c)','Livneh (d, r)','USHCN (m, r)','Berkeley (m, r)'};
dtaStart = [1850 1850 1880 1900 1880 1915 1900 1850];
dtaEnd = [2023 2023 2023 2023 2022 2018 2023 2023];
type = ["m","m","m","m","d","d","m","m"];
var = {'TavgIn','Tin','Tin','TavgIn','TavgIn','TavgIn','TavgIn','TavgIn'};
for d = 1:length(dtasets)
    dtaID = dtasets{d};
    dtaname = dtasetNames{d};
    t = readtable(['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv']);
    y = dtaStart(d)+1:dtaEnd(d);
    T.(dtaname) = t;
    YEARS.(dtaname) = y;
    t.WinterYear = t.Year + (t.Month == 12);
    t.WinterYear(t.Month >= 3 & t.Month <= 11) = NaN;
    if type(d) == "d"
        ssn = seasonFunc(t.(var{d}),@mean,dtaStart(d)+1,dtaEnd(d),dtaStart(d));
        djf = ssn(:,1);
    elseif type(d) == "m"
        djf = nan(length(y),1);
        for yr = 1:length(y)
            djf(yr) = mean(t.(var{d})(t.WinterYear == y(yr)));
        end
    else
        error("Wrong type")
    end

    refI = y >= REF(1) & y <= REF(2);
    djfAnom = djf - mean(djf(refI));
    DJF.(dtaname) = djfAnom;
end




%% Plotting
cols = colororder();
cols(end+1,:) = [0 0 0];
cols = cols([1:3 5 4 8 7 6],:);
figure('Position',[687 533 1000 450]);
hold on
legDta = strings(0,1);
for d = 1:6
    dtaname = dtasetNames{d};
    
    plot(YEARS.(dtaname),DJF.(dtaname),'LineWidth',0.6,'Color',cols(d,:))
    plot(YEARS.(dtaname),runMean(DJF.(dtaname),7),'LineWidth',2,'Color',cols(d,:))
    
    legDta(end+1) = "";
    legDta(end+1) = string(dtasetLegend{d});
end
xlim([1900 2023])
legend(legDta,'location','eastoutside')
set(gca,'FontSize',16,'FontName','Avenir')
xlabel("Year")
ylabel("WH T_{avg} (°C)")
yticks(-6:2:6)
yticklabels(-6:2:6)
grid on
box on
exportgraphics(gcf,'Figures/Supplement/SI1.png','Resolution',450)




hf = figure;
set(hf, 'Position', [50 50 1000 700]);
ha = tight_subplot(2,1,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = [0.08 0.6 0.8 0.35];
ha(2).Position = [0.08 0.1 0.8 0.35];
axes(ha(1))
hold on
legDta = strings(0,1);
for d = [1 4 6 7 8]
    dtaname = dtasetNames{d};
    
    plot(YEARS.(dtaname),DJF.(dtaname),'LineWidth',0.6,'Color',cols(d,:))
    plot(YEARS.(dtaname),runMean(DJF.(dtaname),7),'LineWidth',2,'Color',cols(d,:))
    
    legDta(end+1) = "";
    legDta(end+1) = string(dtasetLegend{d});
end
xlim([1900 2023])
legend(legDta,'location','eastoutside')
set(gca,'FontSize',16,'FontName','Avenir')
xlabel("Year")
ylabel("WH T_{avg} (°C)")
xticks(1900:20:2020)
xticklabels(1900:20:2020)
yticks(-6:2:6)
yticklabels(-6:2:6)
grid on
box on
title("a) Bias-corrected and non-bias-corrected datasets")


axes(ha(2))
hold on
plot(YEARS.BerkeleyMonthly,DJF.BerkeleyMonthly-DJF.BerkeleyRaw,'LineWidth',2,'Color',cols(1,:))
plot(YEARS.USHCNcorr,DJF.USHCNcorr-DJF.USHCNraw,'LineWidth',2,'Color',cols(4,:))
yline(0,'LineWidth',1,'Color','k')

xlim([1900 2023])
legend(["Berkeley","USHCN"],'location','eastoutside')
set(gca,'FontSize',16,'FontName','Avenir')
xlabel("Year")
ylabel("Correction (°C)")
xticks(1900:20:2020)
xticklabels(1900:20:2020)
yticks(-2:0.5:2)
yticklabels(-2:0.5:2)
title("b) Bias correction")
grid on
box on

ha(1).Position = [0.08 0.6 0.75 0.35];
ha(2).Position = [0.08 0.1 0.75 0.35];


exportgraphics(gcf,'Figures/Supplement/SI2.png','Resolution',450)
