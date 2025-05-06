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


function seasons = seasonFuncREP(dta,func,startYear,endYear,seasInds)
    if ~exist('startYear','var')
        % if no startYear is given, assume 1901
        startYear = 1901;
    end
    if ~exist('endYear','var')
        % if no endYear is given, assume 2020
        endYear = 2020;
    end
    
    seasons = zeros(endYear-startYear+1,4);
    
    for year = startYear:endYear
        for ssn = 1:4
            seasons(year-startYear+1,ssn) = func(dta(seasInds(year-startYear+1,1,ssn):seasInds(year-startYear+1,2,ssn)));
        end
    end
end
    