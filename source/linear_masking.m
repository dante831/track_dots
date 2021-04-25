function image = linear_masking(downscale_ratio, theta0, X, Y, image)

    %point1 = [730, 790] * downscale_ratio;
    %point2 = [875, 744] * downscale_ratio;
    
    point1 = downscale_ratio * ([961, 546] + 18 * [cos(theta0/180*pi), sin(theta0/180*pi)]);
    point2 = downscale_ratio * ([961, 546] + 76 * [cos(theta0/180*pi), sin(theta0/180*pi)]);
    

    ind = (Y - point1(2)) * cos((90 + theta0)/180*pi) < (X - point1(1)) * sin((90 + theta0)/180*pi)...
        & (Y - point2(2)) * cos((90 + theta0)/180*pi) > (X - point2(1)) * sin((90 + theta0)/180*pi);
    image(ind) = NaN;
end

