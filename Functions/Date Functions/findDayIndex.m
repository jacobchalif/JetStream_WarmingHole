% Finds the day index

function ind = findDayIndex(year,month,day,startYear)

    if ~exist('startYear','var')
        % if non startYear is given, assume 1901
        startYear = 1901;
    end
    
    yearsCompleted = (startYear:year-1)'; % total number of completed years
    completedLeaps = sum(arrayfun(@(x) isLeap(x),yearsCompleted)); % number of completed leap years
    completedNonLeaps = (year - startYear) - completedLeaps; % number of completed non-leap years
    
    ind = (366*completedLeaps) + (365*completedNonLeaps) + date2yearDay(year,month,day);
end