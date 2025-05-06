WH_latBnds = [30 40];
WH_lonBnds = [-95 -82];
dtaset = 'FLs.52j'; % 'raw' or 'FLs.52j'

%% Read station data
stations = readtable("Data/USHCN/ushcn-v2.5-stations.txt");
NSTATIONS = size(stations,1);
statComplete = true(NSTATIONS,1);
for d = 1:NSTATIONS
    disp(d)
    name = stations.Var1{d};
    [station.(name),completeness] = readTemp(name,dtaset);
    if completeness < 0.8
        station = rmfield(station,name);
        statComplete(d) = false;
        disp(['Removed ' name])
    end
end
stations = stations(statComplete,:);
NSTATIONS = size(stations,1);

%% Interpolation
cell_size = 1;
sigma = 1.5;

LAT = 25:cell_size:49;
LON = -124:cell_size:-67;
YEARS = 1900:2023;
NLAT = length(LAT);
NLON = length(LON);
NYEARS = length(YEARS);
NMONTHS = NYEARS*12;

TMIN = nan(NLAT,NLON,NMONTHS);
TAVG = nan(NLAT,NLON,NMONTHS);
TMAX = nan(NLAT,NLON,NMONTHS);

for lon = 1:NLON
    disp(lon)
    for lat = 1:NLAT
        t = zeros(NMONTHS,3);
        weight_temp = zeros(NMONTHS,3);
        for s = 1:NSTATIONS
            statLat = stations.Var2(s);
            statLon = stations.Var3(s);
            statT = station.(stations.Var1{s});
            latProx = ((LAT(lat) - 1.5*cell_size) <= statLat) && ((LAT(lat) + 1.5*cell_size) >= statLat);
            lonProx = ((LON(lon) - 1.5*cell_size) <= statLon) && ((LON(lon) + 1.5*cell_size) >= statLon);
            if latProx && lonProx
                station_weight = normpdf(0,sigma,sqrt((statLat - LAT(lat))^2 + (statLon - LON(lon))^2));
                statI = find(statT(:,1) >= 1900,1);
                statN = size(statT,1) - statI + 1;
                firstYear = statT(statI,1);
                dateI = find(YEARS == firstYear,1)*12 - 11;
                statT = statT(statI:end,3:5);
                statTnoNaN = statT;
                statTnoNaN(isnan(statTnoNaN)) = 0;
                t(dateI:(dateI+statN-1),:) = t(dateI:(dateI+statN-1),:) + station_weight * statTnoNaN;
                weight_temp(dateI:(dateI+statN-1),:) = weight_temp(dateI:(dateI+statN-1),:) + station_weight * ~isnan(statT);
            end
        end
        if sum(weight_temp,"all") == 0 
            t(:) = NaN;
        end

        TMIN(lat,lon,:) = t(:,1) ./ weight_temp(:,1);
        TAVG(lat,lon,:) = t(:,2) ./ weight_temp(:,2);
        TMAX(lat,lon,:) = t(:,3) ./ weight_temp(:,3);
    end
end

%% Get insideOutside temp
WH_latI = (LAT >= WH_latBnds(1)) & (LAT <= WH_latBnds(2));
WH_lonI = (LON >= WH_lonBnds(1)) & (LON <= WH_lonBnds(2));

tminInside = squeeze(mean(TMIN(WH_latI,WH_lonI,:),[1 2],'omitnan'));
tavgInside = squeeze(mean(TAVG(WH_latI,WH_lonI,:),[1 2],'omitnan'));
tmaxInside = squeeze(mean(TMAX(WH_latI,WH_lonI,:),[1 2],'omitnan'));

tminOutside = squeeze(mean(TMIN(~WH_latI,~WH_lonI,:),[1 2],'omitnan'));
tavgOutside = squeeze(mean(TAVG(~WH_latI,~WH_lonI,:),[1 2],'omitnan'));
tmaxOutside = squeeze(mean(TMAX(~WH_latI,~WH_lonI,:),[1 2],'omitnan'));
dates = nan(NMONTHS,2);
for y = 1:NYEARS
    dates((y-1)*12+[1:12],1) = YEARS(y);
    dates((y-1)*12+[1:12],2) = 1:12;
end
T = [dates tminInside tmaxInside tminOutside tmaxOutside tavgInside tavgOutside];
header = ["Year" "Month" "TminIn" "TmaxIn" "TminOut" "TmaxOut","TavgIn","TavgOut"];
writematrix([header;T],['Data/Processed Reanalysis/SfcAir_InOut/USHCN.' dtaset '.csv'])


function [T,completeness] = readTemp(sID,dtaset)
    minID = fopen(['Data/USHCN/' dtaset '/tmin/' sID '.' dtaset '.tmin']);
    avgID = fopen(['Data/USHCN/' dtaset '/tavg/' sID '.' dtaset '.tavg']);
    maxID = fopen(['Data/USHCN/' dtaset '/tmax/' sID '.' dtaset '.tmax']);
    
    % tmin
    minline = fgetl(minID);
    avgline = fgetl(avgID);
    maxline = fgetl(maxID);
    minyear = str2double(minline(13:16));
    avgyear = str2double(minline(13:16));
    maxyear = str2double(minline(13:16));
    FIRSTYEAR = min([minyear avgyear maxyear]);
    YEARS = FIRSTYEAR:2024;
    NYEARS = length(YEARS);
    T = nan(NYEARS*12,5);
    for y = 1:NYEARS
        T((1:12)+(12*(y-1)),1) = YEARS(y);
        T((1:12)+(12*(y-1)),2) = 1:12;
    end
    mincont = true;
    avgcont = true;
    maxcont = true;
    while mincont || avgcont || maxcont
        if mincont
            minyear = str2double(minline(13:16));
            minI = find(T(:,1) == minyear,1);
            for i = 1:12
                T(minI+i-1,3) = str2double(minline([17:22] + (i-1) * 9));
            end
            minline = fgetl(minID);
            mincont = ischar(minline);
        end

        if avgcont
            avgyear = str2double(avgline(13:16));
            avgI = find(T(:,1) == avgyear,1);
            for i = 1:12
                T(avgI+i-1,4) = str2double(avgline([17:22] + (i-1) * 9));
            end
            avgline = fgetl(avgID);
            avgcont = ischar(avgline);
        end

        if maxcont
            maxyear = str2double(maxline(13:16));
            maxI = find(T(:,1) == maxyear,1);
            for i = 1:12
                T(maxI+i-1,5) = str2double(maxline([17:22] + (i-1) * 9));
            end
            maxline = fgetl(maxID);
            maxcont = ischar(maxline);
        end
    end
    fclose(minID);
    fclose(avgID);
    fclose(maxID);

    LASTYEAR = max([minyear avgyear maxyear]);
    if LASTYEAR == 2024
        LASTYEAR = 2023;
    end
    I = find(T(:,1) == LASTYEAR,1,'last');
    T = T(1:I,:);
    T(T == -9999) = NaN;
    T(:,3:5) = T(:,3:5) / 100;
    if FIRSTYEAR < 1900
        FIRSTYEAR = 1900;
        i = find(T(:,1) == 1900,1);
        T = T(i:end,:);
    end

    %% completeness
    addNaN = 0;
    if FIRSTYEAR > 1900
        addNaN = addNaN + 12*(FIRSTYEAR - 1900);
    end
    if LASTYEAR < 2023
        addNaN = addNaN + 12*(2023 - LASTYEAR);
    end
    completeness = 1 - (sum(sum(isnan(T(:,3:5)),2) > 0) + addNaN) / (124*12);
    
    %% Annual completeness
    for y = FIRSTYEAR:LASTYEAR
        yi = T(:,1) == y;
        t = sum(sum(isnan(T(yi,3:5)),2)>0);
        if t >= 3
            T(yi,3:5) = NaN;
        end
    end
end