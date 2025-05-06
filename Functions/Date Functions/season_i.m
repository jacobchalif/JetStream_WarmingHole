function seasons = season_i(startYR,endYR)
    seasons = false(findDayIndex(endYR,12,31)-findDayIndex(startYR,1,1)+1,4);
    for day = 1:size(seasons,1)
        seasons(day,findSeason(day)) = true;
    end
end