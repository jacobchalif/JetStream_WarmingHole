function yRun = runMean(y,run)
    if ndims(y) > 2
        error("Too many dimensions")
    end
    if size(y,1) == 1
        y = y';
    end
    if size(y,2) > 1
        dim2mean = 2;
    else
        dim2mean = 1;
    end
    step = floor(run/2);
    if (run / 2) == round(run/2)
        error("run must be odd")
    end
    N = length(y);
    yRun = nan(size(y));
    for i = (step+1):(N-step)
        yRun(i,:) = mean(y((i-step):(i+step)),dim2mean);
    end
end