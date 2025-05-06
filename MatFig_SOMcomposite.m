addpath(genpath("Functions"))

SOM_ID = '20240318_NCEP_19802010';
load(['SOMs/' SOM_ID '.mat']);
nclasses = DIM_1*DIM_2;
N = findDayIndex(REF_YEARS(end),12,31,REF_YEARS(1));
dates = nan(N,3);
for year = REF_YEARS
    l = isLeap(year);
    m = monthDays(l);
    i1 = findDayIndex(year,1,1,REF_YEARS(1));
    i2 = findDayIndex(year,12,31,REF_YEARS(1));
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
winterI = (dates(:,2) == 12) | (dates(:,2) <= 2);
dates = dates(winterI,:);

%% Winds
disp("Winds: Reading lat/lon")
winds.lon = ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.2000.nc","lon");
winds.lat = ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.2000.nc","lat");
winds.level = ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.2000.nc","level");
winds.levelI = find(winds.level == LVL);
winds.latI = find(winds.lat >= LAT_MIN & winds.lat <= LAT_MAX);
winds.lonI = find(winds.lon >= LON_MIN & winds.lon <= LON_MAX);
winds.lat = double(winds.lat(winds.latI));
winds.lon = double(winds.lon(winds.lonI))-360;
winds.nlat = length(winds.lat);
winds.nlon = length(winds.lon);
disp("Winds: Reading Winds")
winds.u = nan(winds.nlat,winds.nlon,findDayIndex(REF_YEARS(end),12,31,REF_YEARS(1)));
winds.v = nan(winds.nlat,winds.nlon,findDayIndex(REF_YEARS(end),12,31,REF_YEARS(1)));
for year = REF_YEARS
    i1 = findDayIndex(year,1,1,REF_YEARS(1));
    i2 = findDayIndex(year,12,31,REF_YEARS(1));
    disp([char(num2str(year)) ' u'])
    u = squeeze(ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.' char(num2str(year)) '.nc'],'uwnd',[winds.lonI(1) winds.latI(1) winds.levelI 1],[winds.nlon winds.nlat 1 Inf]));
    winds.u(:,:,i1:i2) = permute(u,[2 1 3]);
    
    disp([char(num2str(year)) ' v'])
    v = squeeze(ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/vwnd_DailyPressure/vwnd.' char(num2str(year)) '.nc'],'vwnd',[winds.lonI(1) winds.latI(1) winds.levelI 1],[winds.nlon winds.nlat 1 Inf]));
    winds.v(:,:,i1:i2) = permute(v,[2 1 3]);
end

winds.u_djf = winds.u(:,:,winterI);
winds.v_djf = winds.v(:,:,winterI);
winds.u_mean = mean(winds.u_djf);
winds.v_mean = mean(winds.v_djf);


%% z
lon = ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.2000.nc","lon");
lat = ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.2000.nc","lat");
lat1 = find(lat==max(lat));
lat2 = find(lat==min(abs(lat)));
level = ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.2000.nc","level");
levelI = find(level == 500); 
latI = find(lat >= LAT_MIN & lat <= LAT_MAX);
lonI = find(lon >= LON_MIN & lon <= LON_MAX);
nlat = length(latI);
nlon = length(lonI);
z = nan(nlon,nlat,findDayIndex(REF_YEARS(end),12,31,REF_YEARS(1)));
for year = REF_YEARS
    i1 = findDayIndex(year,1,1,REF_YEARS(1));
    i2 = findDayIndex(year,12,31,REF_YEARS(1));
    disp([char(num2str(year)) ' LWA'])
    hgtT = squeeze(ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.' char(num2str(year)) '.nc'],'hgt',[1 lat1 levelI 1],[Inf lat2-lat1+1 1 Inf]));
    z(:,:,i1:i2) = hgtT(lonI,latI,:);
end
lat = double(lat(latI));
lon = double(lon(lonI))-360;
z = permute(z,[2 1 3]);
z_djf = z(:,:,winterI);
z_mean = mean(z,3);

%% Finding classes
z_som = reshape(z_djf,[nlat*nlon sum(winterI)]);
y = net(z_som);
classes = vec2ind(y);

z_map = nan(nlat,nlon,nclasses);
u_map = nan(winds.nlat,winds.nlon,nclasses);
v_map = nan(winds.nlat,winds.nlon,nclasses);
for class = 1:nclasses
    i = classes == class;
    z_map(:,:,class) = mean(z_djf(:,:,i),3);
    u_map(:,:,class) = mean(winds.u_djf(:,:,i),3);
    v_map(:,:,class) = mean(winds.v_djf(:,:,i),3);
end

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


%% plotting
baseFS = 12;
titFS = 15;
labFS = 18;

classNames = strings(nclasses,1);
c1 = 0;
c2 = 0;
for class = 1:nclasses
    classNames(class) = strcat("(", num2str(c1), ", " , num2str(c2),")");
    if c2 < (DIM_2-1)
        c2 = c2 + 1;
    else
        c2 = 0;
        c1 = c1 + 1;
    end
end

%states = shaperead('usastatelo', 'UseGeoCoords', true);
land = shaperead('landareas', 'UseGeoCoords', true); %Loads state boundaries


%% Height
figure('Position',[121 343 1100 700]);
ha = tight_subplot(DIM_1,DIM_2,[.05 .01],[.15 .03],[.03 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
CMAX = 5900;
CMIN = 4900;
CSTEP = 100;
for class = 1:nclasses
    axes(ha(class))
    p = contourf(lon,lat,z_map(:,:,class),CMIN:CSTEP:CMAX);
    %p.EdgeColor = 'none';
    hold on
    quiver(lon(1:3:end),lat(1:3:end),u_map(1:3:end,1:3:end,class),v_map(1:3:end,1:3:end,class),'k','LineWidth',1.25);
    clim([CMIN CMAX])
    xlim([LON_MIN LON_MAX]-360)
    ylim([LAT_MIN LAT_MAX])
    if ismember(class,WAVES)
        xl = xlim();
        yl = ylim();
        patch(xl([1 2 2 1 1]),yl([1 1 2 2 1]),'red','FaceColor','none','EdgeColor','red','LineWidth',6)
    end
    %for s = 1:length(states)
    %    plot(states(s).Lon,states(s).Lat,'Color','k','linewidth',1.25)
    %end
    for s = 1:length(land)
        plot(land(s).Lon,land(s).Lat,'Color','k','linewidth',1.25)
    end
    if (class-1)/7 ~= round((class-1)/7)
        yticklabels([])
    else
        yticks(20:10:50)
    end
    if class < 29
        xticklabels([])
    end
    set(gca,'FontSize',baseFS-1,'FontName','Avenir')
    title(classNames(class),'FontSize',titFS)
end
c=contourcbar('southoutside');
c.Position(1)=0.25;
c.Position(2) = 0.08;
c.Position(3) = 0.5;
c.Position(4)= 0.03;
c.Label.String = '500-hPa Geopotential height (m)';
c.FontSize = baseFS;
c.Label.FontSize = labFS;
exportgraphics(gcf,['Figures/Figure3.png'],'Resolution',450)

figure('Position',[100 100 700 50])
c=colormap();
NC = (CMAX-CMIN)/CSTEP;
cols = nan(NC,3);
for i = 1:NC
   x1 = CMIN+(i-1)*CSTEP;
   x2 = CMIN+i*CSTEP;
   cols(i,:) = c(1+floor((i-1)*(length(c)-1)/(NC-1)),:);
   patch([x1 x1 x2 x2 x1],[0 1 1 0 0],cols(i,:))
end
xlim([CMIN CMAX])
xticks(CMIN:CSTEP:CMAX)
yticks([])
set(gca,'FontSize',baseFS,'FontName','Avenir')
exportgraphics(gcf,['Figures/Colorbars/Figure3.png'],'Resolution',450)

%% Winds
figure('Position',[121 343 1100 700]);
ha = tight_subplot(DIM_1,DIM_2,[.05 .01],[.15 .03],[.03 .03]);  % (Nh, Nw, gap, marg_h, marg_w)
CMAX = 40;
CMIN = 0;
CSTEP = 5;
for class = 1:nclasses
    axes(ha(class))
    p = contourf(lon,lat,sqrt(u_map(:,:,class).^2 + v_map(:,:,class).^2),CMIN:CSTEP:CMAX);
    colormap turbo
    %p.EdgeColor = 'none';
    hold on
    quiver(lon(1:3:end),lat(1:3:end),u_map(1:3:end,1:3:end,class),v_map(1:3:end,1:3:end,class),'k','LineWidth',1.25);
    clim([CMIN CMAX])
    xlim([LON_MIN LON_MAX]-360)
    ylim([LAT_MIN LAT_MAX])
    if ismember(class,WAVES)
        xl = xlim();
        yl = ylim();
        patch(xl([1 2 2 1 1]),yl([1 1 2 2 1]),'red','FaceColor','none','EdgeColor','red','LineWidth',6)
    end
    %for s = 1:length(states)
    %    plot(states(s).Lon,states(s).Lat,'Color','k','linewidth',1.25)
    %end
    for s = 1:length(land)
        plot(land(s).Lon,land(s).Lat,'Color','k','linewidth',1.25)
    end
    if (class-1)/7 ~= round((class-1)/7)
        yticklabels([])
    else
        yticks(20:10:50)
    end
    if class < 29
        xticklabels([])
    end
    set(gca,'FontSize',baseFS-1,'FontName','Avenir')
    title(classNames(class),'FontSize',titFS)
end
c=colorbar('southoutside');
c.Position(1)=0.25;
c.Position(2) = 0.08;
c.Position(3) = 0.5;
c.Position(4)= 0.03;
c.Label.String = '300-hPa Wind speed (m/s)';
c.FontSize = baseFS;
c.Label.FontSize = labFS;
exportgraphics(gcf,['Figures/Figure4.png'],'Resolution',450)

figure('Position',[100 100 700 50])
c=colormap('turbo');
NC = (CMAX-CMIN)/CSTEP;
cols = nan(NC,3);
for i = 1:NC
   x1 = CMIN+(i-1)*CSTEP;
   x2 = CMIN+i*CSTEP;
   cols(i,:) = c(1+floor((i-1)*(length(c)-1)/(NC-1)),:);
   patch([x1 x1 x2 x2 x1],[0 1 1 0 0],cols(i,:))
end
xlim([CMIN CMAX])
xticks(CMIN:CSTEP:CMAX)
yticks([])
set(gca,'FontSize',baseFS,'FontName','Avenir')
exportgraphics(gcf,['Figures/Colorbars/Figure4.png'],'Resolution',450)
