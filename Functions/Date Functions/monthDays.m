% Returns vector with number of days in a given month

function months = monthDays(LY)
    if LY
        months = [31,29,31,30,31,30,31,31,30,31,30,31]';  % leap year days
    else
        months = [31,28,31,30,31,30,31,31,30,31,30,31]';  % non-leap year days
    end
end