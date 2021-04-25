function [x0, y0] = find_center(X, Y, x0, y0, r, image)


    % one default solution can be found at 
    % https://www.mathworks.com/matlabcentral/answers/68721-how-to-find-center-point-coordinates-of-circles-in-an-image-file
    % but I won't use that
    [Gmag, ~] = imgradient(image);
    % pick out the most circular part in the plot
    ind_1 = (X - x0).^2 + (Y - y0).^2 > (r*0.9)^2 | (X - x0).^2 + (Y - y0).^2 < (r*0.85)^2;
    ind_2 = Gmag < 25 | Gmag > 150;
    Gmag(ind_1) = NaN;
    Gmag(ind_2) = NaN;
    %figure
    %contour(Gmag1)
    data_Gmag1 = [X(~isnan(Gmag)), Y(~isnan(Gmag))];
    %scatter(data_Gmag1(:, 1), data_Gmag1(:, 2))
    n = 2000;
    center = deal(zeros(n, 2));
    for i = 1 : n
        dopick = true;
        while dopick
            pick_i = ceil(rand(4, 1)*length(data_Gmag1(:, 1)));
            xx = double(data_Gmag1(pick_i, 1));
            yy = double(data_Gmag1(pick_i, 2));
            if (xx(2) - xx(1))^2 + (yy(2) - yy(1))^2 < (r/2)^2 || ...
               (xx(4) - xx(3))^2 + (yy(4) - yy(3))^2 < (r/2)^2
                continue
            end
            sol = intersection(xx, yy);
            if (sol - [x0; y0])' * (sol - [x0; y0]) < (r/5)^2
                center(i, :) = sol;
                dopick = false;
            end
        end
    end
    %scatter(center(:, 1), center(:, 2))
    x0 = mean(center(:, 1));
    y0 = mean(center(:, 2));

    %[IDX, D] = knnsearch(data_Gmag1, data_Gmag1, 'k', 20);
    %data_Gmag1_updated = zeros(size(data_Gmag1));
    %for i = 1 : length(data_Gmag1(:, 1))
    %    data_Gmag1_updated(i, :) = mean(data_Gmag1(IDX(i, :), :));
    %end
    %scatter(data_Gmag1_updated(:, 1), data_Gmag1_updated(:, 2))
    
    
    