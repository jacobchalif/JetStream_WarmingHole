%% Datasets Included: JRA55, NCEP, ERA20c, 20CRv3, ERA5, Livneh

%% Universal Parameters
readLatBnds = [24 52];
readLonBnds = [-125 -65];
WH_latBnds = [30 40];
WH_lonBnds = [-95 -82];
addpath(genpath("Functions"))
setup_nctoolbox
USA_shapefile = 'Data/USAmask/contiguous/contiguous.shp';
USA = shaperead(USA_shapefile);

%% JRA55
dtaID = 'JRA55';
startYear = 1958;
endYear = 2022;
refYear = 1958;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = readLatBnds;
lonBnds = 360+readLonBnds;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/temp_sfc/2021.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/temp_sfc/2021.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgInside = nan(findDayIndex(endYear,12,31,startYear),1);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/temp_sfc/' year '.nc'],"time");
    tsfc=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/temp_sfc/' year '.nc'],"Temperature_height_above_ground",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
    tsfc = permute(squeeze(tsfc) - 273.15,[2 1 3]);
    
    insideT = tsfc;
    insideT(repmat(~WHmask,1,1,size(tsfc,3))) = NaN;
    insideT = squeeze(mean(insideT,[1 2],'omitnan'));
    outsideT = tsfc;
    outsideT(repmat(~outsideMask,1,1,size(tsfc,3))) = NaN;
    outsideT = squeeze(mean(outsideT,[1 2],'omitnan'));
    for t = 1:4:(length(time)-3)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);
        tminInside(ind) = min(insideT(t:t+3));
        tmaxInside(ind) = max(insideT(t:t+3));
        tavgInside(ind) = mean(insideT(t:t+3));
        tminOutside(ind) = min(outsideT(t:t+3));
        tmaxOutside(ind) = max(outsideT(t:t+3));
        tavgOutside(ind) = mean(outsideT(t:t+3));
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
end
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '_no2023.csv'])


%% JRA55 - 2023
tminInside23 = nan(365,1);
tmaxInside23 = nan(365,1);
tavgInside23 = nan(365,1);
tminOutside23 = nan(365,1);
tmaxOutside23 = nan(365,1);
tavgOutside23 = nan(365,1);
dates = nan(365,3);
lastYMD = "";
surf = ncgeodataset('/Volumes/ClimDrive2/Reanalysis/JRA55/2023/surface/anl_surf125.2023010100');
lat=surf.variable('lat'); lat = lat(:);
lon=surf.variable('lon'); lon = lon(:);

latBnds = readLatBnds;
lonBnds = 360+readLonBnds;
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;


for i = 1:365
    [y,m,d] = ind2Date(i,2023);
    YMD = char(strcat("2023",pad(num2str(m),2,'left','0'),pad(num2str(d),2,'left','0')));
    disp(['surface : ' YMD]);
    insideT = nan(4,1);
    outsideT = nan(4,1);
    for h = 1:4
        hourID = char(pad(num2str((h-1)*6),2,'left','0'));
        fileID = ['anl_surf125.' YMD hourID];
        surf = ncgeodataset(['/Volumes/ClimDrive2/Reanalysis/JRA55/2023/surface/' fileID]);
        temp = surf.geovariable('Temperature_height_above_ground');
        temp = squeeze(temp(:,:,:,:) - 273.15);
        insideT(h) = mean(temp(WHmask),'all');
        outsideT(h) = mean(temp(outsideMask),'all');
    end
    tminInside23(i) = min(insideT);
    tmaxInside23(i) = max(insideT);
    tavgInside23(i) = mean(insideT);
    tminOutside23(i) = min(outsideT);
    tmaxOutside23(i) = max(outsideT);
    tavgOutside23(i) = mean(outsideT);
    dates(i,1) = y;
    dates(i,2) = m;
    dates(i,3) = d;
end
T = [dates tminInside23 tmaxInside23 tminOutside23 tmaxOutside23 tavgInside23 tavgOutside23];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
Told = readmatrix('Data/Processed Reanalysis/SfcAir_InOut/JRA55_no2023.csv');
Tnew = [Told;T];
writematrix([header;Tnew],'Data/Processed Reanalysis/SfcAir_InOut/JRA55.csv')



%% NCEP/NCAR
dtaID = 'NCEP';
startYear = 1948;
endYear = 2023;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = readLatBnds;
lonBnds = 360+readLonBnds;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/sfcAir_4xDaily/air.sig995.1948.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/sfcAir_4xDaily/air.sig995.1948.nc","lon");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgOutside = nan(findDayIndex(endYear,12,31,startYear),1);
dates = nan(findDayIndex(endYear,12,31,startYear),3);

for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/sfcAir_4xDaily/air.sig995.' year '.nc'],"time");
    tsfc=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/sfcAir_4xDaily/air.sig995.' year '.nc'],"air",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    tsfc = permute(squeeze(tsfc) - 273.15,[2 1 3]);
    insideT = tsfc;
    insideT(repmat(~WHmask,1,1,size(tsfc,3))) = NaN;
    insideT = squeeze(mean(insideT,[1 2],'omitnan'));
    outsideT = tsfc;
    outsideT(repmat(~outsideMask,1,1,size(tsfc,3))) = NaN;
    outsideT = squeeze(mean(outsideT,[1 2],'omitnan'));
    for t = 1:4:(length(time)-3)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);
        tminInside(ind) = min(insideT(t:t+3));
        tmaxInside(ind) = max(insideT(t:t+3));
        tavgInside(ind) = mean(insideT(t:t+3));
        tminOutside(ind) = min(outsideT(t:t+3));
        tmaxOutside(ind) = max(outsideT(t:t+3));
        tavgOutside(ind) = mean(outsideT(t:t+3));
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
end
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])


%% ERA20c
dtaID = 'ERA20c';
startYear = 1900;
endYear = 2010;
refYear = 1900;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = readLatBnds;
lonBnds = 360+readLonBnds;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA20c/temp_2m/1948.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA20c/temp_2m/1948.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgOutside = nan(findDayIndex(endYear,12,31,startYear),1);
dates = nan(findDayIndex(endYear,12,31,startYear),3);

for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/temp_2m/' year '.nc'],"time");
    tsfc=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/temp_2m/' year '.nc'],"2_metre_temperature_surface",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    tsfc = permute(squeeze(tsfc) - 273.15,[2 1 3]);
    insidetsfc = tsfc;
    insidetsfc(repmat(~WHmask,1,1,size(tsfc,3))) = NaN;
    insidetsfc = mean(insidetsfc,[1 2],'omitnan');
    outsidetsfc = tsfc;
    outsidetsfc(repmat(~outsideMask,1,1,size(tsfc,3))) = NaN;
    outsidetsfc = mean(outsidetsfc,[1 2],'omitnan');
    for t = 1:4:(length(time)-3)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);
        insideT = insidetsfc(t:t+3);
        outsideT = outsidetsfc(t:t+3);
        tminInside(ind) = min(insideT);
        tmaxInside(ind) = max(insideT);
        tavgInside(ind) = mean(insideT);
        tminOutside(ind) = min(outsideT);
        tmaxOutside(ind) = max(outsideT);
        tavgOutside(ind) = mean(outsideT);
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
end
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])



%% ERA5
dtaID = 'ERA5';
startYear = 1940;
endYear = 2023;
refYear = 1900;
YEARS = startYear:endYear;
N = length(YEARS);
fileStarts = [1940 1950 1960 1970 1980 1990 2000 2010 2020 2023];
fileEnds = [1949 1959 1969 1979 1989 1999 2009 2019 2022 2023];
Nfiles = length(fileStarts);
latBnds = readLatBnds;
lonBnds = readLonBnds;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/temp_2m/ERA5_2mTemp_1940_1950.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/temp_2m/ERA5_2mTemp_1940_1950.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;


tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgOutside = nan(findDayIndex(endYear,12,31,startYear),1);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:Nfiles
    yrStart = char(num2str(fileStarts(i)));
    yrEnd = char(num2str(fileEnds(i)+1));
    fileID = ['/Volumes/ClimDrive2/Reanalysis/ERA5/temp_2m/ERA5_2mTemp_' yrStart '_' yrEnd '.nc'];
    lastI = 0;
    for yr = fileStarts(i):fileEnds(i)
        disp([dtaID ': ' char(num2str(yr))])
        time=ncread(fileID,"time",[lastI+1],[(365+isLeap(yr))*4]);
        if yr < 2023
            tsfc=ncread(fileID,"t2m",[lonI(1) latI(1) lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 (365+isLeap(yr))*4]);
        else
            tsfc=ncread(fileID,"t2m",[lonI(1) latI(1) 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 (365+isLeap(yr))*4]);
        end
        tsfc = permute(squeeze(tsfc) - 273.15, [2 1 3]);
        insidetsfc = tsfc;
        insidetsfc(repmat(~WHmask,1,1,size(tsfc,3))) = NaN;
        insidetsfc = mean(insidetsfc,[1 2],'omitnan');
        outsidetsfc = tsfc;
        outsidetsfc(repmat(~outsideMask,1,1,size(tsfc,3))) = NaN;
        outsidetsfc = mean(outsidetsfc,[1 2],'omitnan');
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            insideT = insidetsfc(t:t+3);
            outsideT = outsidetsfc(t:t+3);
            tminInside(ind) = min(insideT);
            tmaxInside(ind) = max(insideT);
            tavgInside(ind) = mean(insideT);
            tminOutside(ind) = min(outsideT);
            tmaxOutside(ind) = max(outsideT);
            tavgOutside(ind) = mean(outsideT);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
        lastI = lastI+(365+isLeap(yr))*4;
    end
end
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])



%% 20CRv3
dtaID = '20CRv3';
startYear = 1900;
endYear = 2015;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = readLatBnds;
lonBnds = 360+readLonBnds;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/tmax_dailySurface/tmax.2m.1900.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/tmax_dailySurface/tmax.2m.1900.nc","lon");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;


tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);

dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/tmax_dailySurface/tmax.2m.' year '.nc'],"time");
    TMIN=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/tmin_dailySurface/tmin.2m.' year '.nc'],"tmin",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    TMAX=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/tmax_dailySurface/tmax.2m.' year '.nc'],"tmax",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    TMIN = permute(squeeze(TMIN - 273.15),[2 1 3]);
    TMAX = permute(squeeze(TMAX - 273.15),[2 1 3]);
    insideTMIN = TMIN;
    insideTMAX = TMAX;
    insideTMIN(repmat(~WHmask,1,1,size(TMIN,3))) = NaN;
    insideTMAX(repmat(~WHmask,1,1,size(TMAX,3))) = NaN;
    insideTMIN = mean(insideTMIN,[1 2],'omitnan');
    insideTMAX = mean(insideTMAX,[1 2],'omitnan');

    outsideTMIN = TMIN;
    outsideTMAX = TMAX;
    outsideTMIN(repmat(~outsideMask,1,1,size(TMIN,3))) = NaN;
    outsideTMAX(repmat(~outsideMask,1,1,size(TMAX,3))) = NaN;
    outsideTMIN = mean(outsideTMIN,[1 2],'omitnan');
    outsideTMAX = mean(outsideTMAX,[1 2],'omitnan');
    for t = 1:length(time)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);
        tminInside(ind) = insideTMIN(t);
        tmaxInside(ind) = insideTMAX(t);
        tminOutside(ind) = outsideTMIN(t);
        tmaxOutside(ind) = outsideTMAX(t);
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
end
tavgInside = mean([tminInside tmaxInside],2);
tavgOutside = mean([tminOutside tmaxOutside],2);
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])


%% Livneh
dtaID = 'Livneh';
startYear = 1915;
endYear = 2018;
YEARS = startYear:endYear;
N = length(YEARS);
lat=ncread("/Volumes/ClimDrive2/Reanalysis/Livneh/tmax/tmax.1915.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/Livneh/tmax/tmax.1915.nc","lon");
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;


tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/Livneh/tmax/tmax.' year '.nc'],"time");
    TMIN=ncread(['/Volumes/ClimDrive2/Reanalysis/Livneh/tmin/tmin.' year '.nc'],"tmin");
    TMAX=ncread(['/Volumes/ClimDrive2/Reanalysis/Livneh/tmax/tmax.' year '.nc'],"tmax");
    TMIN = permute(TMIN,[2 1 3]);
    TMAX = permute(TMAX,[2 1 3]);
    TMIN_WH = TMIN;
    TMIN_WH(repmat(~WHmask,1,1,size(TMIN,3))) = NaN;
    TMIN_WH = squeeze(mean(TMIN_WH,[1 2],'omitnan'));
    TMAX_WH = TMAX;
    TMAX_WH(repmat(~WHmask,1,1,size(TMAX,3))) = NaN;
    TMAX_WH = squeeze(mean(TMAX_WH,[1 2],'omitnan'));

    TMIN_OUTSIDE = TMIN;
    TMIN_OUTSIDE(repmat(~outsideMask,1,1,size(TMIN,3))) = NaN;
    TMIN_OUTSIDE = squeeze(mean(TMIN_OUTSIDE,[1 2],'omitnan'));
    TMAX_OUTSIDE = TMAX;
    TMAX_OUTSIDE(repmat(~outsideMask,1,1,size(TMAX,3))) = NaN;
    TMAX_OUTSIDE = squeeze(mean(TMAX_OUTSIDE,[1 2],'omitnan'));
    for t = 1:length(time)
        d = time(t)+1;
        [y,m,dy] = ind2Date(d,startYear);
        tminInside(d) = TMIN_WH(t);
        tmaxInside(d) = TMAX_WH(t);
        tminOutside(d) = TMIN_OUTSIDE(t);
        tmaxOutside(d) = TMAX_OUTSIDE(t);
        dates(d,1)=y;
        dates(d,2)=m;
        dates(d,3)=dy;
    end
end
tavgInside = mean([tminInside tmaxInside],2);
tavgOutside = mean([tminOutside tmaxOutside],2);
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])





%% GISTEMP
dtaID = 'GISTEMP';
refYear = 1800;
time=ncread("/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/gistemp250_GHCNv4.nc","time");
dates = nan(length(time),2);
for t=1:length(time)
    [y,m] = utConvert(time(t)*24,refYear);
    dates(t,1) = y;
    dates(t,2) = m;
end
dateI = dates(:,1) <= 2023;
dates = dates(dateI,:);
startYear = min(dates(:,1));
endYear = max(dates(:,1));
YEARS = startYear:endYear;
N = length(YEARS);

latBnds = readLatBnds;
lonBnds = readLonBnds;
lat=ncread("/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/gistemp250_GHCNv4.nc","lat");
lon=ncread("/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/gistemp250_GHCNv4.nc","lon");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays
WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

T = ncread("/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/gistemp250_GHCNv4.nc","tempanomaly",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
T = permute(T(:,:,dateI),[2 1 3]);
tInside = T;
tInside(repmat(~WHmask,1,1,size(T,3))) = NaN;
tInside = squeeze(mean(tInside,[1 2],'omitnan'));
tOutside = T;
tOutside(repmat(~outsideMask,1,1,size(T,3))) = NaN;
tOutside = squeeze(mean(tOutside,[1 2],'omitnan'));


T = [dates tInside tOutside];
header = ["Year" "Month" "Tin" "Tout"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])





%% CRUTEM5
dtaID = 'CRUTEM5';
f = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/CRUTEM.5.0.2.0.anomalies.nc";
refYear = 1850;
time=ncread(f,"time");
dates = nan(length(time),2);
for t=1:length(time)
    [y,m] = utConvert(time(t)*24,refYear);
    dates(t,1) = y;
    dates(t,2) = m;
end
dateI = dates(:,1) <= 2023;
dates = dates(dateI,:);
startYear = min(dates(:,1));
endYear = max(dates(:,1));
YEARS = startYear:endYear;
N = length(YEARS);

latBnds = readLatBnds;
lonBnds = readLonBnds;
lat=ncread(f,"latitude");
lon=ncread(f,"longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays
WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

T = ncread(f,"tas",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
T = permute(T(:,:,dateI),[2 1 3]);
tInside = T;
tInside(repmat(~WHmask,1,1,size(T,3))) = NaN;
tInside = squeeze(mean(tInside,[1 2],'omitnan'));
tOutside = T;
tOutside(repmat(~outsideMask,1,1,size(T,3))) = NaN;
tOutside = squeeze(mean(tOutside,[1 2],'omitnan'));


T = [dates tInside tOutside];
header = ["Year" "Month" "Tin" "Tout"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])



%% Berkeley
dtaID = 'Berkeley';
favg = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/Berkeley/Complete_TAVG_LatLong1.nc";
fmin = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/Berkeley/Complete_TMIN_LatLong1.nc";
fmax = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/Berkeley/Complete_TMAX_LatLong1.nc";
refYear = 1850;
time=ncread(fmin,"time");
timeavg=ncread(favg,"time");
dates = nan(length(time),2);
dates(:,1) = floor(time);
dates(1,2) = 1;
for t=2:length(time)
    if dates(t,1) == dates(t-1,1)
        dates(t,2) = dates(t-1,2)+1;
    else
        dates(t,2) = 1;
    end
end
dateI = dates(:,1) <= 2023;
dates = dates(dateI,:);
startYear = min(dates(:,1));
endYear = max(dates(:,1));
YEARS = startYear:endYear;
N = length(YEARS);

avgyears = floor(timeavg);
dateAvgI = avgyears <= 2023 & avgyears >= dates(1,1);
dateAvgI(find(dateAvgI,1,'last')) = false;

latBnds = readLatBnds;
lonBnds = readLonBnds;
lat=ncread(fmin,"latitude");
lon=ncread(fmin,"longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays
WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

TMIN = ncread(fmin,"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
TMAX = ncread(fmax,"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
TAVG = ncread(favg,"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
TMIN = permute(TMIN(:,:,dateI),[2 1 3]);
TMAX = permute(TMAX(:,:,dateI),[2 1 3]);
TAVG = permute(TAVG(:,:,dateAvgI),[2 1 3]);
tminInside = TMIN;
tminInside(repmat(~WHmask,1,1,size(TMIN,3))) = NaN;
tminInside = squeeze(mean(tminInside,[1 2],'omitnan'));
tmaxInside = TMAX;
tmaxInside(repmat(~WHmask,1,1,size(TMAX,3))) = NaN;
tmaxInside = squeeze(mean(tmaxInside,[1 2],'omitnan'));
tavgInside = TAVG;
tavgInside(repmat(~WHmask,1,1,size(TAVG,3))) = NaN;
tavgInside = squeeze(mean(tavgInside,[1 2],'omitnan'));
tminOutside = TMIN;
tminOutside(repmat(~outsideMask,1,1,size(TMIN,3))) = NaN;
tminOutside = squeeze(mean(tminOutside,[1 2],'omitnan'));
tmaxOutside = TMAX;
tmaxOutside(repmat(~outsideMask,1,1,size(TMAX,3))) = NaN;
tmaxOutside = squeeze(mean(tmaxOutside,[1 2],'omitnan'));
tavgOutside = TAVG;
tavgOutside(repmat(~outsideMask,1,1,size(TAVG,3))) = NaN;
tavgOutside = squeeze(mean(tavgOutside,[1 2],'omitnan'));


T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])





%% Berkeley - not QC
dtaID = 'BerkeleyNotQC';
favg = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/Berkeley/NotQC/Raw_TAVG_LatLong1.nc";
fmin = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/Berkeley/NotQC/Raw_TMIN_LatLong1.nc";
fmax = "/Users/jacobchalif/Desktop/WarmingHole/FINAL SCRIPTS/Data/Temperature Analysis/Berkeley/NotQC/Raw_TMAX_LatLong1.nc";
refYear = 1850;
time=ncread(fmin,"time");
timeavg=ncread(favg,"time");
dates = nan(length(time),2);
dates(:,1) = floor(time);
dates(1,2) = 1;
for t=2:length(time)
    if dates(t,1) == dates(t-1,1)
        dates(t,2) = dates(t-1,2)+1;
    else
        dates(t,2) = 1;
    end
end
dateI = dates(:,1) <= 2023;
dates = dates(dateI,:);
startYear = min(dates(:,1));
endYear = max(dates(:,1));
YEARS = startYear:endYear;
N = length(YEARS);

avgyears = floor(timeavg);
dateAvgI = avgyears <= 2023 & avgyears >= dates(1,1);
dateAvgI(find(dateAvgI,1,'last')) = false;

latBnds = readLatBnds;
lonBnds = readLonBnds;
lat=ncread(fmin,"latitude");
lon=ncread(fmin,"longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);

%% Limit to contiguous USA
if max(lon) > 0
    X = repmat((lon-360)',Nlat,1);
else
    X = repmat(lon',Nlat,1);
end
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % Check if the coordinates (X(i), Y(i)) fall within the borders of the USA
    in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    
    % Set the corresponding element in the mask
    USA_mask(i) = in_USA;
end
writematrix(USA_mask,['Data/USAmask/' dtaID '.csv'])

% Inside and Outside WH lat/lon arrays
WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;

TMIN = ncread(fmin,"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
TMAX = ncread(fmax,"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
TAVG = ncread(favg,"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
TMIN = permute(TMIN(:,:,dateI),[2 1 3]);
TMAX = permute(TMAX(:,:,dateI),[2 1 3]);
TAVG = permute(TAVG(:,:,dateAvgI),[2 1 3]);
tminInside = TMIN;
tminInside(repmat(~WHmask,1,1,size(TMIN,3))) = NaN;
tminInside = squeeze(mean(tminInside,[1 2],'omitnan'));
tmaxInside = TMAX;
tmaxInside(repmat(~WHmask,1,1,size(TMAX,3))) = NaN;
tmaxInside = squeeze(mean(tmaxInside,[1 2],'omitnan'));
tavgInside = TAVG;
tavgInside(repmat(~WHmask,1,1,size(TAVG,3))) = NaN;
tavgInside = squeeze(mean(tavgInside,[1 2],'omitnan'));
tminOutside = TMIN;
tminOutside(repmat(~outsideMask,1,1,size(TMIN,3))) = NaN;
tminOutside = squeeze(mean(tminOutside,[1 2],'omitnan'));
tmaxOutside = TMAX;
tmaxOutside(repmat(~outsideMask,1,1,size(TMAX,3))) = NaN;
tmaxOutside = squeeze(mean(tmaxOutside,[1 2],'omitnan'));
tavgOutside = TAVG;
tavgOutside(repmat(~outsideMask,1,1,size(TAVG,3))) = NaN;
tavgOutside = squeeze(mean(tavgOutside,[1 2],'omitnan'));


T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])





%% Berkeley - daily
dtaID = 'BerkeleyDaily';
startYear = 1880;
endYear = 2022;
yearFiles = 1880:10:2020;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = readLatBnds;
lonBnds = readLonBnds;
lat=ncread("Data/Temperature Analysis/Berkeley/Daily/Complete_TAVG_Daily_LatLong1_1900.nc","latitude");
lon=ncread("Data/Temperature Analysis/Berkeley/Daily/Complete_TAVG_Daily_LatLong1_1900.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

if latI(2)<latI(1)
    latI = flip(latI);
    latflip = true;
else
    latflip = false;
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
    lonflip = true;
else
    lonflip = false;
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
Nlat = length(lat);
Nlon = length(lon);


X = repmat(lon',Nlat,1);
Y = repmat(lat,1,Nlon);

% Initialize boolean mask
USA_mask = false(Nlat,Nlon);

frac = 0;
% Loop through each grid cell
for i = 1:numel(X)
    if i/(Nlon*Nlat) - frac > 0.01
        frac = frac + 0.01;
        disp([char(num2str(round(frac*100))) '%'])
    end
    % in_USA = inpolygon(X(i), Y(i), [USA.X], [USA.Y]);
    % USA_mask(i) = in_USA;

    xs = (X(i)-0.5):0.1:(X(i)+0.5);
    ys = (Y(i)-0.5):0.1:(Y(i)+0.5);
    XS = repmat(xs,11,1);
    YS = repmat(ys',1,11);
    in_USA = false(numel(XS),1);
    for i2 = 1:numel(XS)
        in_USA(i2) = inpolygon(XS(i2), YS(i2), [USA.X], [USA.Y]);
    end
    if sum(in_USA)/length(in_USA) > 0.5
        USA_mask(i) = true;
    end
end
writematrix(USA_mask,'Data/USAmask/BerkeleyDaily.csv')


% Inside and Outside WH lat/lon arrays


WHmask = USA_mask & (X >= WH_lonBnds(1)) & (X <= WH_lonBnds(2)) & (Y >= WH_latBnds(1)) & (Y <= WH_latBnds(2));
outsideMask = USA_mask & ~WHmask;


tavgInside = nan(findDayIndex(endYear,12,31,startYear),1);
tavgOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tminInside = nan(findDayIndex(endYear,12,31,startYear),1);
tminOutside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxInside = nan(findDayIndex(endYear,12,31,startYear),1);
tmaxOutside = nan(findDayIndex(endYear,12,31,startYear),1);

dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:length(yearFiles)
    yr = yearFiles(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    ind = findDayIndex(yr,1,1,startYear);
    TAVG=ncread(['Data/Temperature Analysis/Berkeley/Daily/Complete_TAVG_Daily_LatLong1_' year '.nc'],"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    TMIN=ncread(['Data/Temperature Analysis/Berkeley/Daily/Complete_TMIN_Daily_LatLong1_' year '.nc'],"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    TMAX=ncread(['Data/Temperature Analysis/Berkeley/Daily/Complete_TMAX_Daily_LatLong1_' year '.nc'],"temperature",[lonI(1) latI(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 Inf]);
    TAVG = permute(squeeze(TAVG - 273.15),[2 1 3]);
    TMIN = permute(squeeze(TMIN - 273.15),[2 1 3]);
    TMAX = permute(squeeze(TMAX - 273.15),[2 1 3]);
    insideTAVG = TAVG;
    insideTAVG(repmat(~WHmask,1,1,size(TAVG,3))) = NaN;
    insideTAVG = mean(insideTAVG,[1 2],'omitnan');
    outsideTAVG = TAVG;
    outsideTAVG(repmat(~outsideMask,1,1,size(TAVG,3))) = NaN;
    outsideTAVG = mean(outsideTAVG,[1 2],'omitnan');

    insideTMIN = TMIN;
    insideTMIN(repmat(~WHmask,1,1,size(TMIN,3))) = NaN;
    insideTMIN = mean(insideTMIN,[1 2],'omitnan');
    outsideTMIN = TMIN;
    outsideTMIN(repmat(~outsideMask,1,1,size(TMIN,3))) = NaN;
    outsideTMIN = mean(outsideTMIN,[1 2],'omitnan');

    insideTMAX = TMAX;
    insideTMAX(repmat(~WHmask,1,1,size(TMAX,3))) = NaN;
    insideTMAX = mean(insideTMAX,[1 2],'omitnan');
    outsideTMAX = TMAX;
    outsideTMAX(repmat(~outsideMask,1,1,size(TMAX,3))) = NaN;
    outsideTMAX = mean(outsideTMAX,[1 2],'omitnan');
    for t = 1:size(TAVG,3)
        tavgInside(ind+t-1) = insideTAVG(t);
        tavgOutside(ind+t-1) = outsideTAVG(t);
        tmaxInside(ind+t-1) = insideTMAX(t);
        tmaxOutside(ind+t-1) = outsideTMAX(t);
        tminInside(ind+t-1) = insideTMIN(t);
        tminOutside(ind+t-1) = outsideTMIN(t);
    end
end
for i = 1:size(dates,1)
        [y,m,d] = ind2Date(i,1880);           
        dates(i,1)=y;
        dates(i,2)=m;
        dates(i,3)=d;
end
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "Day" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/' dtaID '.csv'])
