function winter = winter_i(startYR,endYR,refYR,rangeEnd)

    if ~exist('refYR','var')
        % if non startYear is given, assume 1901
        refYR = startYR;
    end    
    if ~exist('rangeEnd','var')
        % if non startYear is given, assume 1901
        rangeEnd = endYR;
    end    
    
    winter = false(findDayIndex(rangeEnd,12,31,refYR),1);
    for day = findDayIndex(startYR,1,1,refYR):findDayIndex(endYR,12,31,refYR)
        if findSeason(day) == 1
            winter(day) = true;
        end
    end
end