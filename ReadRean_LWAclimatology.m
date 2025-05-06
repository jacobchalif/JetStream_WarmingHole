
%% NCEP/NCAR
dtaID = 'NCEP';
startYear = 1980;
endYear = 2010;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = map_latBnds;
lonBnds = 360+map_lonBnds;
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


latFinal = lat(latI(1):latI(2));
lonFinal = lon(lonI(1):lonI(2));
nLon = length(lonFinal);
nLat = length(latFinal);
lwa = zeros(nLon,nLat);
nDaysTotal = 0;
for i = 1:N
    yr = YEARS(i);
    days = monthDays(isLeap(yr));
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    hgt=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/hgt_DailyPressure/hgt.' year '.nc'],"hgt",[1 lat1 lvl 1],[Inf lat2-lat1+1 1 Inf]);
    hgt = squeeze(hgt);
    if latflip
        hgt = flip(hgt,2);
    end
    LWA = calculateLWA(hgt,lat,lon);
    winterI = [1:sum(days(1:2)) (sum(days(1:11)+1):sum(days))];
    LWA_ant = LWA(lonI(1):lonI(2),latI(1):latI(2),winterI,2);
    LWA_cyc = LWA(lonI(1):lonI(2),latI(1):latI(2),winterI,3);
    
    lwa = lwa + sum(LWA_ant - LWA_cyc,3);
    nDaysTotal = nDaysTotal + length(winterI);
end

lwa = lwa' / nDaysTotal;

writematrix(lwa,['Data/Processed Reanalysis/Gridded/NCEP/LWA/climatology_1980_2010.csv'])
writematrix(latFinal,['Data/Processed Reanalysis/Gridded/NCEP/LWA/lat.csv'])
writematrix(lonFinal,['Data/Processed Reanalysis/Gridded/NCEP/LWA/lon.csv'])



