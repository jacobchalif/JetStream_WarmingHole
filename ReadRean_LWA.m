%% Datasets Included: JRA55, NCEP, ERA20c, 20CRv3, ERA5

%% Universal Parameters
WH_latBnds = [30 50];
WH_lonBnds = [-100 -60];
glob_latBnds = [30 65];
glob_lonBnds = [-180 180];
levels = {'500'};
addpath(genpath("Functions"))
setup_nctoolbox

latBnds2Use = glob_latBnds;
lonBnds2Use = glob_lonBnds;
reg = 'global';

%% NCEP/NCAR
dtaID = 'NCEP';
startYear = 1948;
endYear = 2023;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
level=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.1948.nc","level");
lat=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.1948.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.1948.nc","lon");
lvl = find(level == 500);
latEq = find(lat==min(abs(lat)));
latPole = find(lat==max(lat));
if lat(1)<lat(2)
    latflip = true;
    lat2 = latPole;
    lat1 = latEq;
    lat = flip(lat(lat1:lat2));
    % lat = flip(lat);
else
    latflip = false;
    lat1 = latPole;
    lat2 = latEq;
    lat = lat(lat1:lat2);
end
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end
LWA_ant = nan(length(lon),length(lat),findDayIndex(endYear,12,31,startYear));
LWA_cyc = nan(length(lon),length(lat),findDayIndex(endYear,12,31,startYear));
lwa = nan(findDayIndex(endYear,12,31,startYear),3);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.' year '.nc'],"time");
    hgt=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.' year '.nc'],"hgt",[1 lat1 lvl 1],[Inf lat2-lat1+1 1 Inf]);
    hgt = squeeze(hgt);
    if latflip
        hgt = flip(hgt,2);
    end
    LWA = calculateLWA(hgt,lat,lon);
    for t = 1:length(time)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);
        lwa(ind,:) = mean(LWA(lonI(1):lonI(2),latI(1):latI(2),t,:),[1 2]);
        LWA_ant(:,:,ind) = LWA(:,:,t,2);
        LWA_cyc(:,:,ind) = LWA(:,:,t,3);
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
end
T = [dates lwa(:,2) lwa(:,3)];
header = ["Year" "Month" "Day" "LWA_Anticylonic" "LWA_Cylonic"];
writematrix([header;T],['Data/Processed Reanalysis/LWA_' reg '/' dtaID '.csv'])


%% JRA55
dtaID = 'JRA55';
startYear = 1958;
endYear = 2022;
refYear = 1958;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/500hPa/hgt/2021.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/500hPa/hgt/2021.nc","longitude");
latEq = find(lat==min(abs(lat)));
latPole = find(lat==max(lat));
if lat(1)<lat(2)
    latflip = true;
    lat2 = latPole;
    lat1 = latEq;
    lat = flip(lat(lat1:lat2));
else
    latflip = false;
    lat1 = latPole;
    lat2 = latEq;
    lat = lat(lat1:lat2);
end
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end

dates = nan(findDayIndex(endYear,12,31,startYear),3);
lwa = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/500hPa/hgt/' year '.nc'],"time");
    HGT=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/500hPa/hgt/' year '.nc'],"Geopotential_height_isobaric_surface_low",[1 lat1 1 1],[Inf lat2-lat1+1 1 Inf]);
    HGT = squeeze(HGT);
    hgt = nan(size(HGT,1),size(HGT,2),size(HGT,3)/4);
    for t = 1:4:(length(time)-3)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);   
        hgt(:,:,(t-1)/4+1) = mean(HGT(:,:,t:t+3),3);
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
    LWA = calculateLWA(hgt,lat,lon);
    startInd = findDayIndex(yr,1,1,startYear);
    endInd = findDayIndex(yr,12,31,startYear);
    lwa(startInd:endInd,:) = mean(LWA(lonI(1):lonI(2),latI(1):latI(2),:,:),[1 2]);
end
T = [dates lwa(:,2) lwa(:,3)];
header = ["Year" "Month" "Day" "LWA_Anticylonic" "LWA_Cylonic"];
writematrix([header;T],['Data/Processed Reanalysis/LWA_' reg '/JRA55_no2023.csv'])


%% JRA55 - 2023
s.lat = [0 90];
dates = nan(365,3);
lastYMD = "";
fil = ncgeodataset('/Volumes/ClimDrive2/Reanalysis/JRA55/2023/hgt/anl_p125_hgt.2023010100');
h=fil.geovariable('Geopotential_height_isobaric');
level=fil.geovariable('isobaric');
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lon=fil.geovariable('lon');
lon=lon(:);
lat=fil.geovariable('lat');
lat=lat(:);
lvl = find(level(:) == 500);
if string(reg) == "global"
    [rowmin, rowmax, colmin, colmax] = geoij(h, s);
    colmin = 1;
    si = size(h);
    colmax = si(4);
else
    [rowmin, rowmax, colmin, colmax] = geoij(h, s);
end
lat = lat(rowmin:rowmax);
if lat(1)<lat(2)
    latflip = true;
    lat = flip(lat(lat1:lat2));
else
    latflip = false;
    lat = lat(lat1:lat2);
end
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));
lwa23 = nan(365,3);
for i = 1:365
    [y,m,d] = ind2Date(i,2023);
    YMD = char(strcat("2023",pad(num2str(m),2,'left','0'),pad(num2str(d),2,'left','0')));
    disp(['LWA : ' YMD]);
    hgt = nan(rowmax-rowmin+1,colmax-colmin+1,4);
    for h = 1:4
        hourID = char(pad(num2str((h-1)*6),2,'left','0'));
        fileID = ['anl_p125_hgt.' YMD hourID];
        fil = ncgeodataset(['/Volumes/ClimDrive2/Reanalysis/JRA55/2023/hgt/' fileID]);
        hgtT = fil.geovariable('Geopotential_height_isobaric');
        hgt(:,:,h) = squeeze(hgtT(1,lvl,rowmin:rowmax,colmin:colmax));
    end
    hgt = mean(hgt,3);
    hgt = hgt';
    if latflip
        hgt = flip(hgt,2);
    end
    LWA = squeeze(calculateLWA(hgt,lat,lon));
    lwa23(i,:) = mean(LWA(lonI(1):lonI(2),latI(1):latI(2),:),[1 2]);
    dates(i,1) = y;
    dates(i,2) = m;
    dates(i,3) = d;
end
T = [dates lwa23(:,2) lwa23(:,3)];
header = ["Year" "Month" "Day" "LWA_Anticylonic" "LWA_Cylonic"];
Told = readmatrix(['Data/Processed Reanalysis/LWA_' reg '/JRA55_no2023.csv']);
Tnew = [Told;T];
writematrix([header;Tnew],['Data/Processed Reanalysis/LWA_' reg '/JRA55.csv'])


%% 20CRv3
dtaID = '20CRv3';
startYear = 1900;
endYear = 2015;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
level=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/hgt_dailyPressure/hgt.1979.nc","level");
lat=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/hgt_dailyPressure/hgt.1979.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/hgt_dailyPressure/hgt.1979.nc","lon");
lvl = find(level == 500);
latEq = find(lat==min(abs(lat)));
latPole = find(lat==max(lat));
if lat(1)<lat(2)
    latflip = true;
    lat2 = latPole;
    lat1 = latEq;
    lat = flip(lat(lat1:lat2));
    % lat = flip(lat);
else
    latflip = false;
    lat1 = latPole;
    lat2 = latEq;
    lat = lat(lat1:lat2);
end
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end
lwa = nan(findDayIndex(endYear,12,31,startYear),3);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/hgt_dailyPressure/hgt.' year '.nc'],"time");
    hgt=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/hgt_dailyPressure/hgt.' year '.nc'],"hgt",[1 lat1 lvl 1],[Inf lat2-lat1+1 1 Inf]);
    hgt = squeeze(hgt);
    if latflip
        hgt = flip(hgt,2);
    end
    LWA = calculateLWA(hgt,lat,lon);
    for t = 1:length(time)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);
        lwa(ind,:) = mean(LWA(lonI(1):lonI(2),latI(1):latI(2),t,:),[1 2]);
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
end
T = [dates lwa(:,2) lwa(:,3)];
header = ["Year" "Month" "Day" "LWA_Anticylonic" "LWA_Cylonic"];
writematrix([header;T],['Data/Processed Reanalysis/LWA_' reg '/' dtaID '.csv'])


%% ERA20c
dtaID = 'ERA20c';
startYear = 1900;
endYear = 2010;
refYear = 1900;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA20c/500hPa/hgt/1900.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA20c/500hPa/hgt/1900.nc","longitude");
latEq = find(lat==min(abs(lat)));
latPole = find(lat==max(lat));
if lat(1)<lat(2)
    latflip = true;
    lat2 = latPole;
    lat1 = latEq;
    lat = flip(lat(lat1:lat2));
    % lat = flip(lat);
else
    latflip = false;
    lat1 = latPole;
    lat2 = latEq;
    lat = lat(lat1:lat2);
end
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end
lwa = nan(findDayIndex(endYear,12,31,startYear),3);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:N
    yr = YEARS(i);
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    time=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/500hPa/hgt/' year '.nc'],"time");
    geo=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/500hPa/hgt/' year '.nc'],"Geopotential_isobaric",[1 lat1 1 1],[Inf lat2-lat1+1 1 Inf]);
    geo = squeeze(geo);
    if latflip
        geo = flip(geo,2);
    end
    hgt = nan(size(geo,1),size(geo,2),size(geo,3)/4);
    for t = 1:4:(length(time)-3)
        d = time(t);
        [y,m,dy,hrs] = utConvert(d,refYear);
        ind = findDayIndex(y,m,dy,startYear);   
        hgt(:,:,(t-1)/4+1) = mean(geo(:,:,t:t+3)/9.80665,3);
        dates(ind,1)=y;
        dates(ind,2)=m;
        dates(ind,3)=dy;
    end
    LWA = calculateLWA(hgt,lat,lon);
    startInd = findDayIndex(yr,1,1,startYear);
    endInd = findDayIndex(yr,12,31,startYear);
    lwa(startInd:endInd,:) = mean(LWA(lonI(1):lonI(2),latI(1):latI(2),:,:),[1 2]);
end
T = [dates lwa(:,2) lwa(:,3)];
header = ["Year" "Month" "Day" "LWA_Anticylonic" "LWA_Cylonic"];
writematrix([header;T],['Data/Processed Reanalysis/LWA_' reg '/' dtaID '.csv'])



%% ERA5
dtaID = 'ERA5';
startYear = 1940;
endYear = 2023;
refYear = 1900;
stride = 2;
YEARS = startYear:endYear;
N = length(YEARS);
fileStarts = [1940 1950 1960 1970 1980 1990 2000 2010 2020 2023];
fileEnds = [1949 1959 1969 1979 1989 1999 2009 2019 2022 2023];
Nfiles = length(fileStarts);
latBnds = latBnds2Use;
lonBnds = lonBnds2Use;
level=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/hgt/ERA5_hgt_300_500_hPa_1940_1950.nc","level");
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/hgt/ERA5_hgt_300_500_hPa_1940_1950.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/hgt/ERA5_hgt_300_500_hPa_1940_1950.nc","longitude");
lvl = find(level == 500);
latEq = find(lat==min(abs(lat)));
latPole = find(lat==max(lat));
if lat(1)<lat(2)
    latflip = true;
    lat2 = latPole;
    lat1 = latEq;
    lat = flip(lat(lat1:lat2));
    % lat = flip(lat);
else
    latflip = false;
    lat1 = latPole;
    lat2 = latEq;
    lat = lat(lat1:lat2);
end
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));

stride = 4;
lon = lon(1:stride:end-(stride-1));
lat = lat(1:stride:end);
[~,latI(1)] = min(abs(lat-latBnds(2)));
[~,latI(2)] = min(abs(lat-latBnds(1)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end

lwa = nan(findDayIndex(endYear,12,31,startYear),3);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:Nfiles
    yrStart = char(num2str(fileStarts(i)));
    yrEnd = char(num2str(fileEnds(i)+1));
    file = ['/Volumes/ClimDrive2/Reanalysis/ERA5/hgt/ERA5_hgt_300_500_hPa_' yrStart '_' yrEnd '.nc'];
    lastI = 0;
    for yr = fileStarts(i):fileEnds(i)
        time=ncread(file,"time",[lastI+1],[(365+isLeap(yr))*4]);
        disp([dtaID ': ' char(num2str(yr))])
        if yr < 2023
            geo=ncread(file,"z",[1 lat1 lvl lastI+1],[Inf (lat2-lat1)/stride+1 1 (365+isLeap(yr))*4],[stride stride 1 1]);
        else
            geo=ncread(file,"z",[1 lat1 lvl 1 lastI+1],[Inf (lat2-lat1)/stride+1 1 1 (365+isLeap(yr))*4],[stride stride 1 1 1]);
        end     
        
        geo = squeeze(geo);
        if latflip
            geo = flip(geo,2);
        end
        hgt = nan(size(geo,1),size(geo,2),size(geo,3)/4);
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);   
            hgt(:,:,(t-1)/4+1) = mean(geo(:,:,t:t+3)/9.80665,3);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
        LWA = calculateLWA(hgt,lat,lon);
        startInd = findDayIndex(yr,1,1,startYear);
        endInd = findDayIndex(yr,12,31,startYear);
        lwa(startInd:endInd,:) = mean(LWA(lonI(1):lonI(2),latI(1):latI(2),:,:),[1 2]);
        lastI = lastI+(365+isLeap(yr))*4;
    end
end
T = [dates lwa(:,2) lwa(:,3)];
header = ["Year" "Month" "Day" "LWA_Anticylonic" "LWA_Cylonic"];
writematrix([header;T],['Data/Processed Reanalysis/LWA_' reg '/' dtaID '.csv'])
