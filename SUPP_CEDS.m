ceds = readmatrix("Data/CEDS_Emissions/SO2.csv");
years = ceds(:,1);
so2 = ceds(:,2);

figure('Position',[584 581 560 275])
plot(years,so2/1e9,'color','k','LineWidth',2)
xlabel("Year")
ylabel("SO_2 emissions (Tg)")
xlim([years(1) years(end)])

set(gca,'FontSize',16,'FontName','Avenir')
exportgraphics(gcf,'Figures/Supplement/SI6.png','Resolution',450)


