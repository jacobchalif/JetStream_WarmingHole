function [year,month,day,hours] = utConvert(hour,startYear)
    hour = double(hour);
    ttime = hour/24;
    days_elap = floor(ttime) + 1;
    hours = 24*(ttime - floor(ttime));
    
    [year,month,day] = ind2Date(days_elap,startYear);
    
end