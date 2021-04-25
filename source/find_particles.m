function particles = find_particles(image, x0_m, y0_m, ub, lb, r_particle)
    
    %ub = 80;
    %lb = 20;
    %r_particle = 3;
    
    ind = image < ub & image > lb;
    
    particles = [];
    
    [y_, x_] = find(ind);
    x_reserve = x_;
    y_reserve = y_;
    
    while ~isempty(x_)
        %particles(length(particles) + 1).pos = [x_(1), y_(1)];
%        if isempty(x_(2:end))
%            particles(length(particles) + 1).points = [x_(1), y_(1)];
%            x_(1) = [];
%            y_(1) = [];
        findpoint = true;
        particles(length(particles) + 1).x_ = x_(1);
        particles(end).y_ = y_(1);
        while findpoint == true && ~isempty(x_)
            temp_x = particles(end).x_;
            temp_y = particles(end).y_;
            dis = (x_' - temp_x) .* (x_' - temp_x) + (y_' - temp_y) .* (y_' - temp_y);
            indices = dis < 3;
            new_points = find(sum(double(indices), 1) >= 1);
            if isempty(new_points)
                findpoint = false;
            else
                particles(end).x_ = [particles(end).x_; x_(new_points)];
                particles(end).y_ = [particles(end).y_; y_(new_points)];
                x_(new_points) = [];
                y_(new_points) = [];
            end
        end
    end
    
    % filter the particles
    i = 1;
    while i <= length(particles)
        temp_x = particles(i).x_;
        temp_y = particles(i).y_;
        dis = (temp_x' - temp_x) .* (temp_x' - temp_x) + (temp_y' - temp_y) .* (temp_y' - temp_y);
        if max(dis(:)) > 10 * r_particle^2 || length(temp_x) < 6
            particles(i) = [];
        else
            particles(i).pos = [mean(particles(i).x_); mean(particles(i).y_)];
            i = i + 1;
        end
    end
    
    %figure
    %scatter(x_reserve - x0_m, y_reserve - y0_m)
    %hold on
    %for i = 1 : length(particles)
        %scatter(particles(i).x_ - x0_m, particles(i).y_ - y0_m)
    %end
    %hold off

    
end
        
        
        
