lwa = readmatrix("Data/Processed Reanalysis/Gridded/NCEP/LWA/climatology_1980_2010.csv");
mci = readmatrix("Data/Processed Reanalysis/Gridded/NCEP/MCI/climatology_1980_2010.csv");
u = readmatrix("Data/Processed Reanalysis/Gridded/NCEP/U_climatology_1980_2010.csv");
v = readmatrix("Data/Processed Reanalysis/Gridded/NCEP/V_climatology_1980_2010.csv");
lat = readmatrix("Data/Processed Reanalysis/Gridded/NCEP/MCI/lat.csv");
lon = readmatrix("Data/Processed Reanalysis/Gridded/NCEP/MCI/lon.csv");
lon = lon - 360;


cellsize=2.5;
skip = 3;

%% Creates figure
hf = figure;
set(hf, 'Position', [50 50 1000 450]);
ha = tight_subplot(1,2,[.04 .05],[.06 .06],[.04 .04]);  % (Nh, Nw, gap, marg_h, marg_w)

axes(ha(1))

h=contourf(lon,flip(lat),flip(lwa)/1e7);
xlim([-150 -30]);
ylim([10 77.5]);
colormap(redblue)
A = colorbar('southoutside'); % put colorbar at the bottom
A.Label.String = strcat("LWA (10^7 m^2)"); % name colorbar title
A.FontSize = 12; 
A.Label.FontSize = 14; 
rectangle('Position',[-100 30 40 20],'LineWidth',2.5,'LineStyle','--','EdgeColor','#ff001e')
hold on 
quiver(lon(1:skip:end),lat(1:skip:end),u(1:skip:end,1:skip:end),v(1:skip:end,1:skip:end),'k','LineWidth',1);
        
caxis([-5 5])
grid off
shading interp
load coastlines
geoshow(coastlat,coastlon,'Color','k')
states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(states,'FaceColor',[1,1,1],'facealpha',0); %Adds state boundaries 
hold on 
    
set(gca,'XTick',-180:30:0,'XTickLabel',-180:30:0)
set(gca,'YTick',0:15:90,'YTickLabel',0:15:90)
title("a: LWA Winter Climatology, 1980-2010",'fontsize',15)
set(gca,'FontName','Avenir')


axes(ha(2))
h=contourf(lon,flip(lat),flip(mci));
xlim([-150 -30]);
ylim([10 77.5]);
colormap(redblue)
A = colorbar('southoutside'); % put colorbar at the bottom
A.Label.String = strcat("MCI"); % name colorbar title
A.FontSize = 12; 
A.Label.FontSize = 14; 
rectangle('Position',[-110 30 20 20],'LineWidth',2.5,'LineStyle','--','EdgeColor','#ff001e')
hold on
quiver(lon(1:skip:end),lat(1:skip:end),u(1:skip:end,1:skip:end),v(1:skip:end,1:skip:end),'k','LineWidth',1);
        
caxis([-0.5,0.5])
grid off
shading interp
load coastlines
geoshow(coastlat,coastlon,'Color','k')
states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(states,'FaceColor',[1,1,1],'facealpha',0); %Adds state boundaries 
hold on 
set(gca,'XTick',-180:30:0,'XTickLabel',-180:30:0)
set(gca,'YTick',0:15:90,'YTickLabel',0:15:90)
set(gca,'FontName','Avenir')
title("b: MCI Winter Climatology, 1980-2010",'fontsize',15)


exportgraphics(gcf,['Figures/Figure1.png'],'Resolution',450)


