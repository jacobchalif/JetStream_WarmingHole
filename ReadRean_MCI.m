%% Datasets Included: JRA55, NCEP, ERA20c, 20CRv3, ERA5

%% Universal Parameters
WH_latBnds = [30 50];
WH_lonBnds = [-110 -90];
glob_latBnds = [30 65];
glob_lonBnds = [0 360];
levels = {'300','500'};
addpath(genpath("Functions"))
setup_nctoolbox

latBnds2Use = glob_latBnds;
lonBnds2Use = glob_lonBnds;
reg = 'global';

%% JRA55
dtaID = 'JRA55';
startYear = 1958;
endYear = 2022;
refYear = 1958;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = latBnds2Use;
lonBnds = 360+lonBnds2Use;
lat=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/300hPa/vwnd/2021.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/JRA55/300hPa/vwnd/2021.nc","longitude");
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
    mci = nan(findDayIndex(endYear,12,31,startYear),1);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ': ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/' levels{l} 'hPa/vwnd/' year '.nc'],"time");
        v=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/' levels{l} 'hPa/vwnd/' year '.nc'],"v-component_of_wind_isobaric_surface_low",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/JRA55/' levels{l} 'hPa/uwnd/' year '.nc'],"u-component_of_wind_isobaric_surface_low",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        MCI = v .* abs(v) ./ (u.^2 + v.^2);
        if string(reg) == "global"
            MCI = squeeze(mean(abs(MCI),[1 2 3]));
        else
            MCI = squeeze(mean(MCI),[1 2 3]);
        end
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            M = MCI(t:t+3);
            mci(ind) = mean(M);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates mci];
    header = ["Year" "Month" "Day" "MCI"];
    writematrix([header;T],['Data/Processed Reanalysis/MCI' levels{l} '_' reg '/' dtaID '_no2023.csv'])
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
    [rowmin, rowmax, colmin, colmax] = geoij(u, s);
    colmin = 1;
    si = size(u);
    colmax = si(3);
else
    [rowmin, rowmax, colmin, colmax] = geoij(u, s);
end
mci23 = nan(365,2);
for i = 1:365
    [y,m,d] = ind2Date(i,2023);
    YMD = char(strcat("2023",pad(num2str(m),2,'left','0'),pad(num2str(d),2,'left','0')));
    disp(['MCI ' levels{l} ' : ' YMD]);
    mciT = nan(4,2);
    for h = 1:4
        hourID = char(pad(num2str((h-1)*6),2,'left','0'));
        ufileID = ['anl_p125_ugrd.' YMD hourID];
        vfileID = ['anl_p125_vgrd.' YMD hourID];
        Ufil = ncgeodataset(['/Volumes/ClimDrive2/Reanalysis/JRA55/2023/ugrd/' ufileID]);
        u = Ufil.geovariable('u-component_of_wind_isobaric');
        Vfil = ncgeodataset(['/Volumes/ClimDrive2/Reanalysis/JRA55/2023/vgrd/' vfileID]);
        v = Vfil.geovariable('v-component_of_wind_isobaric');
        uT = squeeze(u(1,lvls,rowmin:rowmax,colmin:colmax));
        vT = squeeze(v(1,lvls,rowmin:rowmax,colmin:colmax));
        if string(reg) == "global"
            mciT(h,:) = mean(abs(vT .* abs(vT) ./ (uT.^2 + vT.^2)),[2 3]);
        else
            mciT(h,:) = mean(vT .* abs(vT) ./ (uT.^2 + vT.^2),[2 3]);
        end
    end
    mci23(i,:) = mean(mciT,1);
    dates(i,1) = y;
    dates(i,2) = m;
    dates(i,3) = d;
end
for l = 1:2
    T = [dates mci23(:,l)];
    header = ["Year" "Month" "Day" "MCI"];
    Told = readmatrix(['Data/Processed Reanalysis/MCI' levels{l}  '_' reg '/JRA55_no2023.csv']);
    Tnew = [Told;T];
    writematrix([header;Tnew],['Data/Processed Reanalysis/MCI' levels{l} '_' reg '/JRA55.csv'])
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
    mci = nan(findDayIndex(endYear,12,31,startYear),1);
    dates = nan(findDayIndex(endYear,12,31,startYear),3);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ', ' levels{l} ' hPa: ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/vwnd_DailyPressure/vwnd.' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.' year '.nc'],"uwnd",[lonI(1) latI(1) lvls(l) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        v=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/vwnd_DailyPressure/vwnd.' year '.nc'],"vwnd",[lonI(1) latI(1) lvls(l) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        if string(reg) == "global"
            MCI = squeeze(mean(abs(v .* abs(v) ./ (u.^2+v.^2)),[1 2]));
        else
            MCI = squeeze(mean(v .* abs(v) ./ (u.^2+v.^2),[1 2]));
        end
        
        for t = 1:length(time)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            mci(ind) = MCI(t);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates mci];
    header = ["Year" "Month" "Day" "MCI"];
    writematrix([header;T],['Data/Processed Reanalysis/MCI' levels{l} '_' reg '/' dtaID '.csv'])
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
    mci = nan(findDayIndex(endYear,12,31,startYear),1);
    dates = nan(findDayIndex(endYear,12,31,startYear),3);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ': ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/' levels{l} 'hPa/uwnd/' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/' levels{l} 'hPa/uwnd/' year '.nc'],"U_component_of_wind_isobaric",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        v=ncread(['/Volumes/ClimDrive2/Reanalysis/ERA20c/' levels{l} 'hPa/vwnd/' year '.nc'],"V_component_of_wind_isobaric",[lonI(1) latI(1) 1 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        if string(reg) == "global"
            MCI = squeeze(mean(abs(v .* abs(v) ./ (u.^2+v.^2)),[1 2]));
        else
            MCI = squeeze(mean(v .* abs(v) ./ (u.^2+v.^2),[1 2]));
        end
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            M = MCI(t:t+3);
            mci(ind) = mean(M);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates mci];
    header = ["Year" "Month" "Day" "MCI"];
    writematrix([header;T],['Data/Processed Reanalysis/MCI' levels{l} '_' reg '/' dtaID '.csv'])
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
lat=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/vwnd/ERA5_vwnd_300_500_hPa_1940_1950.nc","latitude");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/ERA5/vwnd/ERA5_vwnd_300_500_hPa_1940_1950.nc","longitude");
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
mci = nan(findDayIndex(endYear,12,31,startYear),2);
dates = nan(findDayIndex(endYear,12,31,startYear),3);
for i = 1:Nfiles
    yrStart = char(num2str(fileStarts(i)));
    yrEnd = char(num2str(fileEnds(i)+1));
    vFile = ['/Volumes/ClimDrive2/Reanalysis/ERA5/vwnd/ERA5_vwnd_300_500_hPa_' yrStart '_' yrEnd '.nc'];
    uFile = ['/Volumes/ClimDrive2/Reanalysis/ERA5/uwnd/ERA5_uwnd_300_500_hPa_' yrStart '_' yrEnd '.nc'];
    lastI = 0;
    for yr = fileStarts(i):fileEnds(i)
        time=ncread(vFile,"time",[lastI+1],[(365+isLeap(yr))*4]);
        disp([dtaID ', v: ' char(num2str(yr))])
        if yr < 2023
            v=ncread(vFile,"v",[lonI(1) latI(1) 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 2 (365+isLeap(yr))*4]);
        else
            v=ncread(vFile,"v",[lonI(1) latI(1) 1 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 2 1 (365+isLeap(yr))*4]);
        end
        disp([dtaID ', u: ' char(num2str(yr))])
        if yr < 2023
            u=ncread(uFile,"u",[lonI(1) latI(1) 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 2 (365+isLeap(yr))*4]);
        else
            u=ncread(uFile,"u",[lonI(1) latI(1) 1 1 lastI+1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 2 1 (365+isLeap(yr))*4]);
        end        
        if string(reg) == "global"
            MCI = squeeze(mean(abs(v .* abs(v) ./ (u.^2+v.^2)),[1 2]));
        else
            MCI = squeeze(mean(v .* abs(v) ./ (u.^2+v.^2),[1 2]));
        end
        MCI = MCI';
        for t = 1:4:(length(time)-3)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            M = MCI(t:t+3,:);
            mci(ind,:) = mean(M,1);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
        lastI = lastI+(365+isLeap(yr))*4;
    end
end
for l = 1:2
    T = [dates mci(:,l)];
    header = ["Year" "Month" "Day" "MCI"];
    writematrix([header;T],['Data/Processed Reanalysis/MCI' levels{l} '_' reg '/' dtaID '.csv'])
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
    mci = nan(findDayIndex(endYear,12,31,startYear),1);
    dates = nan(findDayIndex(endYear,12,31,startYear),3);
    for i = 1:N
        yr = YEARS(i);
        year = char(num2str(yr));
        disp([dtaID ', ' levels{l} ': ' year])
        time=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.' year '.nc'],"time");
        u=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/uwnd_dailyPressure/uwnd.' year '.nc'],"uwnd",[lonI(1) latI(1) lvls(l) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        v=ncread(['/Volumes/ClimDrive2/Reanalysis/20CRv3/vwnd_dailyPressure/vwnd.' year '.nc'],"vwnd",[lonI(1) latI(1) lvls(l) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
        if string(reg) == "global"
            MCI = squeeze(mean(abs(v .* abs(v) ./ (u.^2+v.^2)),[1 2]));
        else
            MCI = squeeze(mean(v .* abs(v) ./ (u.^2+v.^2),[1 2]));
        end
        for t = 1:length(time)
            d = time(t);
            [y,m,dy,hrs] = utConvert(d,refYear);
            ind = findDayIndex(y,m,dy,startYear);
            mci(ind) = MCI(t);
            dates(ind,1)=y;
            dates(ind,2)=m;
            dates(ind,3)=dy;
        end
    end
    T = [dates mci];
    header = ["Year" "Month" "Day" "MCI"];
    writematrix([header;T],['Data/Processed Reanalysis/MCI' levels{l} '_' reg '/' dtaID '.csv'])
end
