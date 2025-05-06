% Checks if a year 'year' is a leap year
% Returns true in a variable 'leap' if it is
% Returns false in a variable 'leap' if it is not

function leap = isLeap(year)
    if mod(year,4) ~= 0
        leap = false;
    elseif mod(year,100) ~= 0
        leap = true;
    elseif mod(year,400) ~= 0
        leap = false;
    else
        leap = true;
    end
end

