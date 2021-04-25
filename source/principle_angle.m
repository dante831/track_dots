function theta = principle_angle(X, Y, x0_m, y0_m, r, image_shifted, image0_shifted, string)

    image1 = image_shifted;
    image0 = image0_shifted;
    
    [Gmag, Gdir] = imgradient(image1);
    % mask out wherever is outside r or inside r/3.5
    ind_1 = (X - x0_m).^2 + (Y - y0_m).^2 > (r*0.79)^2 | (X - x0_m).^2 + (Y - y0_m).^2 < (r/3.0)^2;
    %ind_2 = Gmag < 25 | Gmag > 200;
    ind_2 = Gmag < 25 | Gmag > 150;
    ind_3 = image1 < 100;
    Gdir(ind_1) = NaN;
    Gdir(ind_2) = NaN;
    image1(ind_3) = NaN;
    
    [Gmag0, Gdir0] = imgradient(image0);
    ind_2_0 = Gmag0 < 25 | Gmag0 > 150;
    ind_3_0 = image0 < 100;
    Gdir(ind_1) = NaN;
    Gdir(ind_2_0) = NaN;
    image0(ind_3_0) = NaN;
    
    if strcmp(string, 'sample')
        
        % sampling

        error = zeros(360, 1);
        gradient0 = imgradient(image0);
        gradient1 = imgradient(image1);
        ind0 = gradient0 < 100 | gradient0 > 300;
        gradient0(ind0) = 0;
        for i = 1 : length(error)
            %imgradient(image0)>100 .* image0>5.0*mean(mean(image0))            
            gradient1_rotated = imrotate(gradient1, i, 'crop');
            %ind1 = double(gradient1_rotated > 100 & gradient1_rotated < 300);
            ind1 = gradient1_rotated < 100 | gradient1_rotated > 300;
            gradient1_rotated(ind1) = 0;
            error(i) = (gradient0(:) - gradient1_rotated(:))' * (gradient0(:) - gradient1_rotated(:));
            %error(i) = sum(sum(abs(double(image0>5.5*mean(mean(image0))) - double(imrotate(image, i, 'crop')>5.5*mean(mean(image1))))));
            %error(i) = sum(sum(double(abs(image0 - imrotate(image, i, 'crop')))>40));
        end
        %figure
        %plot(error);
        
        angle = (-1:1/100:1) + find(error == min(error));
        error2 = zeros(size(angle));
        for i = 1 : length(angle)
            gradient1_rotated = imrotate(gradient1,  angle(i), 'crop');
            %ind1 = double(gradient1_rotated > 100 & gradient1_rotated < 300);
            ind1 = gradient1_rotated < 100 | gradient1_rotated > 300;
            gradient1_rotated(ind1) = 0;
            error2(i) = (gradient0(:) - gradient1_rotated(:))' * (gradient0(:) - gradient1_rotated(:));
            %error2(i) = (ind0(:) - ind1(:))' * (ind0(:) - ind1(:));
            %error(i) = sum(sum(abs(double(imgradient(image0)>100 & imgradient(image0)<300) - double(imrotate(imgradient(image1), angle(i), 'crop')>100 & imrotate(imgradient(image1), i, 'crop') < 300))));
            %error(i) = sum(sum(double(abs(imgradient(image0)>5.5*mean(mean(image0)) - imrotate(image, i, 'crop')>5.5*mean(mean(image))))));
            %error2(i) = sum(sum(double(abs(image0 - imrotate(image, angle(i), 'crop')))>40));
        end
        %figure
        %plot(error2);
        
        theta = mean(angle(error2 == min(error2)));
        
        
    elseif strcmp(string, 'linear_regression')
        
        % linear regression
        
        temp_ind = ~isnan(Gdir);
        y_new = Y(temp_ind(:));
        x_new = [X(temp_ind(:)), ones(size(y_new))];
        A = [x_new, y_new];
        b = regress(y_new, x_new);
        theta = acos(b(1))/pi*180;
        
        % remove outliers and do regression again to reduce bias of the
        % black dots
        ind_3 = (y_new - x_new * b).^2 > (0.1*r)^2;
        y_new_2 = y_new(ind_3);
        x_new_2 = x_new(ind_3, :);
        %b = regress(y_new_2, x_new_2);
        %theta = acos(b(1))/pi*180;
        
        figure;
        hold on
        scatter(x_new(:, 1), y_new)
        scatter(x_new(:, 1), x_new * b)
        
    elseif strcmp(string, 'PCA')
        
        % PCA
    
        temp_ind = ~isnan(Gdir);
        x_new = X(temp_ind(:));
        y_new = Y(temp_ind(:));
        A = [x_new, y_new];
        [U,S,V] = svd(A);
        theta = acos(V(1, 1))/pi*180;
        %theta2 = acos(V(2, 1))/pi*180;
    
    elseif strcmp(string, 'gradient')
    
        % gradient

        %Gdir((sum(double(Gdir(:) ==  [-180, -90, 0, 90]), 2) == 1)) = NaN;
        [N, edges] = histcounts(Gdir(:), 7200);
        [~, N_ind] = sort(N, 'descend');
        % two directions of the bar
        if abs(abs((edges(N_ind(1)) + edges(N_ind(1) + 1)) / 2 - (edges(N_ind(2)) + edges(N_ind(2) + 1)) / 2) - 180) < 0.1
            theta = ((edges(N_ind(1)) + edges(N_ind(1) + 1)) / 2 + ...
                     (edges(N_ind(2)) + edges(N_ind(2) + 1)) / 2 + 180) / 2;
        else
            theta = (edges(N_ind(1)) + edges(N_ind(1) + 1)) / 2;
        end

        figure;
        contour(Gdir)
        figure;
        polarhistogram(Gdir(:)/180*pi, 7200);
    
    end

