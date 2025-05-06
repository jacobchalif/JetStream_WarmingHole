% Finds the season a particular day is in

function season = findSeason(ind,startYear)
    if ~exist('startYear','var')
        % if non startYear is given, assume 1901
        startYear = 1901;
    end
    
    [~,month,~] = ind2Date(ind,startYear);
    
    if (month <= 2) || (month == 12)
        season = 1;
    elseif month <= 5
        season = 2;
    elseif month <= 8
        season = 3;
    else
        season = 4;
    end
end
        