%% Datasets Included: JRA55, NCEP, ERA20c, 20CRv3, ERA5

%% Universal Parameters
WH_latBnds = [40 60];
WH_lonBnds = [-140 0];
glob_latBnds = [30 65];
glob_lonBnds = [0 360];
levels = {'300','500'};
addpath(genpath("Functions"))
setup_nctoolbox

latBnds2Use = WH_latBnds;
lonBnds2Use = WH_lonBnds;
reg = 'WH';

%% JRA55
dtaID = 'JRA55';
startYear = 1958;
endYear = 2022;
refYear = 1958;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/300hPa/uwnd/2021.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/300hPa/uwnd/2021.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end
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
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for l = 1:2
    uwnd = nan(findDayIndex(endYear,12,31,startYear),1);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ': ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/' levels{l} 'hPa/uwnd/' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/' levels{l} 'hPa/uwnd/' year '.nc'],"u-component_of_wind_isobaric_surface_low",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        u = squeeze(mean(u,[1 2 3]));
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            U = u(t:t+3);
            uwnd(ind) = mean(U);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates uwnd];
    header = ["Year" "Month" "Day" "UWnd"];
    if string(reg) == "global"
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '_' reg '/' dtaID '_no2023.csv'])
    else
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '/' dtaID '_no2023.csv'])
    end
end

%% JRA55 - 2023
s.lat = latBnds2Use;
s.lon = lonBnds2Use;
dates = nan(365,3);
lastYMD = "";
fil = ncgeodataset('/Volumes/ClimDrive2/Reanalysis/JRA55/2023/ugrd/anl_p125_ugrd.2023010100');
u=fil.geovariable('u-component_of_wind_isobaric');
level=fil.geovariable('isobaric');
lvls = [find(level(:) == 300) find(level(:) == 500)];
if string(reg) == "global"
    [rowmin, rowmax, ~, ~] = geoij(u, s);
    colmin = 1;
    si = size(u);
    colmax = si(3);
else
    [rowmin, rowmax, ~, ~] = geoij(u, s);
    dtaLon = u.getlondata();
    [~,colmax] = min(abs(dtaLon-s.lon(1)));
    [~,colmin] = min(abs(dtaLon-s.lon(2)));
end
u23 = nan(365,2);
for i = 1:365
    [y,m,d] = ind2Date(i,2023);
    YMD = char(strcat("2023",pad(num2str(m),2,'left','0'),pad(num2str(d),2,'left','0')));
    disp(['U ' levels{l} ' : ' YMD]);
    uwndT = nan(4,2);
    for h = 1:4
        hourID = char(pad(num2str((h-1)*6),2,'left','0'));
        ufileID = ['anl_p125_ugrd.' YMD hourID];
        Ufil = ncgeodataset(['/Volumes/ClimDrive2/Reanalysis/JRA55/2023/ugrd/' ufileID]);
        u = Ufil.geovariable('u-component_of_wind_isobaric');
        uT = squeeze(u(1,lvls,rowmin:rowmax,colmin:colmax));
        uwndT(h,:) = mean(uT,[2 3]);
    end 
    u23(i,:) = mean(uwndT,1);
    dates(i,1) = y;
    dates(i,2) = m;
    dates(i,3) = d;
end
for l = 1:2
    T = [dates u23(:,l)];
    header = ["Year" "Month" "Day" "UWnd"];
    if string(reg) == "global"
        Told = readmatrix(['Data/Processed Reanalysis/U' levels{l}  '_' reg '/JRA55_no2023.csv']);
    else
        Told = readmatrix(['Data/Processed Reanalysis/U' levels{l} '/JRA55_no2023.csv']);
    end
    Tnew = [Told;T];
    if string(reg) == "global"
        writematrix([header;Tnew],['Data/Processed Reanalysis/U' levels{l} '_' reg '/JRA55.csv'])
    else
        writematrix([header;Tnew],['Data/Processed Reanalysis/U' levels{l} '/JRA55.csv'])
    end
end


%% NCEP/NCAR
dtaID = 'NCEP';
startYear = 1948;
endYear = 2023;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
level=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.1948.nc","level");
lat=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.1948.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.1948.nc","lon");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end
lvls = [find(level == 300) find(level == 500)];

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
for l = 1:2
    uwnd = nan(findDayIndex(endYear,12,31,startYear),1);
    dates = nan(findDayIndex(endYear,12,31,startYear),3);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ', ' levels{l} ' hPa: ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/vwnd_DailyPressure/vwnd.' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.' year '.nc'],"uwnd",[lonI(1) latI(1) lvls(l) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        u = squeeze(mean(u,[1 2]));
        
        for t = 1:length(time)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            uwnd(ind) = u(t);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates uwnd];
    header = ["Year" "Month" "Day" "UWnd"];
    if string(reg) == "global"
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '_' reg '/' dtaID '.csv'])
    else
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '/' dtaID '.csv'])
    end
end


%% ERA20c
dtaID = 'ERA20c';
startYear = 1900;
endYear = 2010;
refYear = 1900;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA20c/500hPa/uwnd/1900.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA20c/500hPa/uwnd/1900.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end

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
for l = 1:2
    uwnd = nan(findDayIndex(endYear,12,31,startYear),1);
    dates = nan(findDayIndex(endYear,12,31,startYear),3);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ': ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/' levels{l} 'hPa/uwnd/' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/' levels{l} 'hPa/uwnd/' year '.nc'],"U_component_of_wind_isobaric",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        u = squeeze(mean(u,[1 2]));
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            U = u(t:t+3);
            uwnd(ind) = mean(U);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates uwnd];
    header = ["Year" "Month" "Day" "UWnd"];
    if string(reg) == "global"
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '_' reg '/' dtaID '.csv'])
    else
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '/' dtaID '.csv'])
    end
end


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
latBnds = latBnds2Use;
lonBnds = lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/uwnd/ERA5_uwnd_300_500_hPa_1940_1950.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/uwnd/ERA5_uwnd_300_500_hPa_1940_1950.nc","longitude");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end

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
uwnd = nan(findDayIndex(endYear,12,31,startYear),2);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:Nfiles
    yrStart = char(num2str(fileStarts(i)));
    yrEnd = char(num2str(fileEnds(i)+1));
    uFile = ['/Volumes/ClimDrive2/Reanalysis/ERA5/uwnd/ERA5_uwnd_300_500_hPa_' yrStart '_' yrEnd '.nc'];
    lastI = 0;
    for yr = fileStarts(i):fileEnds(i)
        time=ncread(uFile,"time",[lastI+1],[(365+isLeap(yr))*4]);
        disp([dtaID ', u: ' char(num2str(yr))])
        if yr < 2023
            u=ncread(uFile,"u",[lonI(1) latI(1) 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 2 (365+isLeap(yr))*4]);
        else
            u=ncread(uFile,"u",[lonI(1) latI(1) 1 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 2 1 (365+isLeap(yr))*4]);
        end        
        u = squeeze(mean(u,[1 2]));
        u = u';
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            U = u(t:t+3,:);
            uwnd(ind,:) = mean(U,1);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
        lastI = lastI+(365+isLeap(yr))*4;
    end
end
for l = 1:2
    T = [dates uwnd(:,l)];
    header = ["Year" "Month" "Day" "UWnd"];
    if string(reg) == "global"
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '_' reg '/' dtaID '.csv'])
    else
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '/' dtaID '.csv'])
    end
end


%% 20CRv3
dtaID = '20CRv3';
startYear = 1900;
endYear = 2015;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.1900.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.1900.nc","lon");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
if string(reg) == "global"
    [~,lonI(1)] = min(lon);
    [~,lonI(2)] = max(lon);
else
    [~,lonI(1)] = min(abs(lon-lonBnds(1)));
    [~,lonI(2)] = min(abs(lon-lonBnds(2)));
end
level=ncread("/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.1900.nc","level");
lvls = [find(level == 300) find(level == 500)];

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
for l = 1:2
    uwnd = nan(findDayIndex(endYear,12,31,startYear),1);
    dates = nan(findDayIndex(endYear,12,31,startYear),3);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ', ' levels{l} ': ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.' year '.nc'],"uwnd",[lonI(1) latI(1) lvls(l) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        u = squeeze(mean(u,[1 2]));
        for t = 1:length(time)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            uwnd(ind) = u(t);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates uwnd];
    header = ["Year" "Month" "Day" "UWnd"];
    if string(reg) == "global"
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '_' reg '/' dtaID '.csv'])
    else
        writematrix([header;T],['Data/Processed Reanalysis/U' levels{l} '/' dtaID '.csv'])
    end
end
