function seasons = seasonInds(startYear,endYear)
    if ~exist('startYear','var')
        % if no startYear is given, assume 1901
        startYear = 1901;
    end
    if ~exist('endYear','var')
        % if no endYear is given, assume 2020
        endYear = 2020;
    end
    
    seasons = zeros(endYear-startYear+1,2,4);
    
    for year = startYear:endYear
        if year == startYear
            seasons(1,1,1) = 1;
            seasons(1,2,1) = findDayIndex(year,3,1,startYear)-1;
        else
            seasons(year-startYear+1,1,1) = findDayIndex(year-1,12,1,startYear);
            seasons(year-startYear+1,2,1) = findDayIndex(year,3,1,startYear)-1;
        end
        
        seasons(year-startYear+1,1,2) = findDayIndex(year,3,1,startYear);
        seasons(year-startYear+1,2,2) = findDayIndex(year,5,31,startYear);
        seasons(year-startYear+1,1,3) = findDayIndex(year,6,1,startYear);
        seasons(year-startYear+1,2,3) = findDayIndex(year,8,31,startYear);
        seasons(year-startYear+1,1,4) = findDayIndex(year,9,1,startYear);
        seasons(year-startYear+1,2,4) = findDayIndex(year,11,30,startYear);
    end
end

