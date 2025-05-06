% does reverse of findDayIndex

function [year,month,day] = ind2Date(ind,startYear)
    if ~exist('startYear','var')
        % if non startYear is given, assume 1901
        startYear = 1901;
    end
    
    day = ind;
    month = 1;
    year = startYear;
    isYear = false;
    isMonth = false;
    
    while ~isYear
        leap = isLeap(year);
        if ~leap && (day > 365)
            year = year + 1;
            day = day - 365;
        elseif leap && (day > 366)
            year = year + 1;
            day = day - 366;
        else
            isYear = true;
        end
    end
    
    months = monthDays(isLeap(year));
    while ~isMonth
        if day > months(month)
            day = day - months(month);
            month = month + 1;
        else
            isMonth = true;
        end
    end
    
end
            
        