% Function to find seasonal values after being passed some function 'func'
    % (returns nyears x 4 array)
% to pass function, include @ before function
    % for example:
        % seasonMeans = seasonFunc(dta,@mean)
        % OR
        % seasonSTD = seasonFunc(dta,@std)
% winter = DJF (first column)
% spring = MAM (second column)
% summer = JJA (third column)
% fall = SON (fourth column)


function seasons = seasonFunc(dta,func,startYear,endYear,refYear)
    if ~exist('startYear','var')
        % if no startYear is given, assume 1901
        startYear = 1901;
    end
    if ~exist('endYear','var')
        % if no endYear is given, assume 2020
        endYear = 2020;
    end
    if ~exist('refYear','var')
        % if no endYear is given, assume 2020
        refYear = startYear;
    end
    seasons = zeros(endYear-startYear+1,4);
    
    for year = startYear:endYear
        if (year == startYear) && (refYear == startYear)
            seasons(1,1) = func(dta(1:findDayIndex(year,3,1,refYear)-1));
        else
            seasons(year-startYear+1,1) = func(dta(findDayIndex(year-1,12,1,refYear):findDayIndex(year,3,1,refYear)-1));
        end
        
        seasons(year-startYear+1,2) = func(dta(findDayIndex(year,3,1,refYear):findDayIndex(year,5,31,refYear)));
        seasons(year-startYear+1,3) = func(dta(findDayIndex(year,6,1,refYear):findDayIndex(year,8,31,refYear)));
        seasons(year-startYear+1,4) = func(dta(findDayIndex(year,9,1,refYear):findDayIndex(year,11,30,refYear)));
    end
end
    