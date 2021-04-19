function [dist, iCount, jCount] = circleSpacingFromNumber(traces, angle, number)
    TraceLimits = getRotatedMapLimits(traces, angle);

    x00 = TraceLimits(1,1);
    x01 = TraceLimits(1,2);
    x10 = TraceLimits(1,3);

    y00 = TraceLimits(2,1);
    y01 = TraceLimits(2,2);
    y10 = TraceLimits(2,3);

    Li = hypot((x01 - x00),(y01 - y00));
    Lj = hypot((x10 - x00),(y10 - y00));
    if Li >= Lj        
        jCount = number;
        dist = Lj / (jCount + 1);
        iCount = floor((Li / dist) - 1);        
    %% Y dimension is longer than x
    else       
        iCount = number;
        dist = Li / (iCount + 1);
        jCount = floor((Lj / dist) - 1);        
    end
end

