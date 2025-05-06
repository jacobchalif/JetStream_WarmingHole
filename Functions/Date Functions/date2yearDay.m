% Function to return index of date in a given year

function ind = date2yearDay(year,month,day)
    days = monthDays(isLeap(year)); % stores number of days in each month
    ind = day; % adds days passed in current month to index;
    ind = ind + sum(days(1:month-1)); % adds number of days that have passed before current month
end
