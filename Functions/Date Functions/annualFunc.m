% Function to find annual values after being passed some function 'func'
% to pass function, include @ before function
    % for example:
        % annualMeans = annualFunc(dta,@mean)
        % OR
        % annualSTD = annualFunc(dta,@std)
% returns nyears x 1 vector

function annual = annualFunc(dta,func,startYear,endYear)
    if ~exist('startYear','var')
        % if no startYear is given, assume 1901
        startYear = 1901;
    end
    if ~exist('endYear','var')
        % if no endYear is given, assume 2020
        endYear = 2020;
    end
    
    annual = zeros(endYear - startYear + 1,1);
    
    for year = startYear:endYear
        yearInds = findDayIndex(year,1,1,startYear):findDayIndex(year,12,31,startYear);
        annual(year-startYear+1,1) = func(dta(yearInds));
    end
end