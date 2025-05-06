addpath(genpath("Functions/"))


mov = 7;
lw = 1.25;
YL = [-3.5 3.5];
titFS = 18;
legFS = 14;
labFS = 17;
baseFS = 12;
        

dtasets = {'NCEP','20CRv3','ERA5','ERA20c','JRA55'};
dtaLeg = {'NCEPr1','20CRv3','ERA5','ERA20c','JRA55'};
legDta = dtaLeg;
starts = [1948 1900 1940 1900 1958];
ends = [2023 2015 2023 2010 2023];
cols = colororder();
cols(end+1,:) = [0 0 0];
overlap = [max(starts) min(ends)];
overlap = [1980 2010];
ssn = 1; % 1 = DJF




%% just timeseries
hf = figure;
set(hf, 'Position', [50 50 500 600]);
ha = tight_subplot(2,1,[.1 .08],[.135 .06],[.06 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
ha(1).Position = ([0.15 0.55 0.72 0.4]);
ha(2).Position = ([0.15 0.1 0.72 0.4]);


for d = 1:length(dtasets)
    disp(dtasets{d})
    years = starts(d)+1:ends(d);
    yearOverlap = years>=overlap(1) & years<=overlap(2);

    % panel b: LWA
    lwa=readtable(['Data/Processed Reanalysis/LWA/' dtasets{d} '.csv']);
    lwa_ssn = seasonFunc(lwa.LWA_Anticylonic+lwa.LWA_Cylonic,@mean,starts(d)+1,ends(d),starts(d));
    lwa_anom = lwa_ssn - mean(lwa_ssn(yearOverlap,:),1);

    axes(ha(1))
    plot(years,runMean(lwa_anom(:,ssn)/1e6,mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir','YDir','Reverse')

    % panel c: MCI
    mci=readtable(['Data/Processed Reanalysis/MCI300_global/' dtasets{d} '.csv']);
    mci_ssn = seasonFunc(mci.MCI,@mean,starts(d)+1,ends(d),starts(d));
    mci_anom = mci_ssn - mean(mci_ssn(yearOverlap,:),1);

    axes(ha(2))
    plot(years,runMean(mci_anom(:,ssn),mov),'LineWidth',lw,'Color',cols(d,:))
    hold on
    set(gca,'FontSize',baseFS,'FontName','Avenir','YDir','Reverse')
    xlabel("Year",'FontSize',labFS)
end
for i = 1:2
    axes(ha(i))
    xticks(1900:20:2020)
    xlim([1900 2023])
    grid on
    if i < 2
        xticklabels([])
    end
end

axes(ha(1))
ylabel("DJF LWA_{NH} anomaly (10^{6} m^{2})",'FontSize',labFS)
text(2020,-3.4,"a: LWA — Northern Hemisphere",'FontSize',titFS,'FontName','Avenir','HorizontalAlignment','right')
ylim([-4 4])
yticks(-4:2:4)

axes(ha(2))
ylabel("DJF MCI_{NH} anomaly",'FontSize',labFS)
text(2020,-0.014,"b: MCI — Northern Hemisphere",'FontSize',titFS,'FontName','Avenir','HorizontalAlignment','right')
ylim([-0.016 0.016])
yticks(-0.01:0.01:0.01)
legend(legDta,'NumColumns',4,'Location','southwest','FontSize',legFS)

exportgraphics(gcf,['Figures/Figure8.png'],'Resolution',450)
