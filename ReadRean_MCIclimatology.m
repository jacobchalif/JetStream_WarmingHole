%% Universal Parameters
WH_latBnds = [30 50];
WH_lonBnds = [-110 -90];
map_latBnds = [15 75];
map_lonBnds = [-150 -30];
addpath(genpath("Functions"))

%% NCEP/NCAR
dtaID = 'NCEP';
startYear = 1980;
endYear = 2010;
refYear = 1800;
YEARS = startYear:endYear;
N = length(YEARS);
latBnds = map_latBnds;
lonBnds = 360+map_lonBnds;
level=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.1948.nc","level");
lat=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.1948.nc","lat");
lon=ncread("/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.1948.nc","lon");
[~,latI(1)] = min(abs(lat-latBnds(1)));
[~,latI(2)] = min(abs(lat-latBnds(2)));
[~,lonI(1)] = min(abs(lon-lonBnds(1)));
[~,lonI(2)] = min(abs(lon-lonBnds(2)));
lvls = find(level == 300);

if latI(2)<latI(1)
    latI = flip(latI);
end
if lonI(2)<lonI(1)
    lonI = flip(lonI);
end
lat = lat(latI(1):latI(2));
lon = lon(lonI(1):lonI(2));
nLon = length(lon);
nLat = length(lat);
mci = zeros(nLon,nLat);
u = zeros(nLon,nLat);
v = zeros(nLon,nLat);
nDaysTotal = 0;
for i = 1:N
    yr = YEARS(i);
    days = monthDays(isLeap(yr));
    year = char(num2str(yr));
    disp([dtaID ': ' year])
    U=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/uwnd_DailyPressure/uwnd.' year '.nc'],"uwnd",[lonI(1) latI(1) lvls(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
    V=ncread(['/Volumes/ClimDrive2/Reanalysis/NCEP/vwnd_DailyPressure/vwnd.' year '.nc'],"vwnd",[lonI(1) latI(1) lvls(1) 1],[lonI(2)-lonI(1)+1 latI(2)-latI(1)+1 1 Inf]);
    winterI = [1:sum(days(1:2)) (sum(days(1:11)+1):sum(days))]; 
    U = U(:,:,winterI);
    V = V(:,:,winterI);
    MCI = squeeze(V .* abs(V) ./ (U.^2+V.^2));
    mci = mci + sum(MCI,3);
    u = u + sum(U,3);
    v = v + sum(V,3);
    nDaysTotal = nDaysTotal + length(winterI);
end
mci = mci' / nDaysTotal;
u = u' / nDaysTotal;
v = v' / nDaysTotal;

writematrix(mci,['Data/Processed Reanalysis/Gridded/NCEP/MCI/climatology_1980_2010.csv'])
writematrix(u,['Data/Processed Reanalysis/Gridded/NCEP/U_climatology_1980_2010.csv'])
writematrix(v,['Data/Processed Reanalysis/Gridded/NCEP/V_climatology_1980_2010.csv'])
writematrix(lat,['Data/Processed Reanalysis/Gridded/NCEP/MCI/lat.csv'])
writematrix(lon,['Data/Processed Reanalysis/Gridded/NCEP/MCI/lon.csv'])




