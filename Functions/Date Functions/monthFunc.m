% Function to find monthly values after being passed some function 'func'
    % to pass function, include @ before function
    % for example:
        % monthMeans = monthFunc(dta,@mean)
        % OR
        % monthSTD = monthfunc(dta,@std)
% if grid = false (default), returns 12*nyears x 3 vector
    % first column = year
    % second column = month
    % third column = monthly value
% if grid = true, returns nyears x 12 array

function months = monthFunc(dta,func,grid,startYear,endYear)
    if ~exist('grid','var')
        % if no grid is given, assume sequential timeseries
        grid = false;
    end
    if ~exist('startYear','var')
        % if non startYear is given, assume 1901
        startYear = 1901;
    end
    if ~exist('endYear','var')
        % if no endYear is given, assume 2020
        endYear = 2020;
    end
    
    if grid
        months = zeros(startYear-endYear+1,12);
    else
        months = zeros((startYear-endYear+1)*12,3);
    end
    
    for year = startYear:endYear
        for month = 1:12
            monthInds = findDayIndex(year,month,1):findDayIndex(year,month+1,1)-1;
            if grid
                months(year-startYear+1,month) = func(dta(monthInds));
            else
                months((year-startYear)*12+month,1) = year;
                months((year-startYear)*12+month,2) = month;
                months((year-startYear)*12+month,3) = func(dta(monthInds));
            end
        end
    end
end