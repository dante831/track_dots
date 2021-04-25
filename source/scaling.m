

    % select precipitations from 07021800-07041730
    
    precip_file_start = 1;
    precip_file_end = 432;
    precip_files = precip_files(precip_file_start : precip_file_end);

    % select lat-lon region

    lat_series = h5read(precip_files(1).name, '/Grid/lat');
    lon_series = h5read(precip_files(1).name, '/Grid/lon');
    dphi = lat_series(2) - lat_series(1);
    dlambda = lon_series(2) - lon_series(1);
    
    [lat_indices, lon_indices] = latlonindices(lat_series, lon_series, latmin, latmax, lonmin, lonmax);   
    lat = lat_series(lat_indices);
    lon = lon_series(lon_indices);
    
    % defining and tracking algorithms
    % first, find the precipitation event based on threshold
    %precip_events = struct('start_time', [], 'end_time', [], ...
    %                       'region_lon', [], 'region_lat', [], ...
    %                       'mean_lon', [], 'mean_lat', [], ...
    %                       'added', [], 'accu_precip', []);  
                            % 'added' is a logical parameter, being 1 if this event is added
                            %       with points. This parameter is used to determine whether 
                            %       this event has ended.
    precip_events = [];
    precip_events_output = [];
    
    lat_series_2 = ncread(nc_filename, 'latitude');
    lat_series_2 = lat_series_2(end:-1:1);
    lon_series_2 = ncread(nc_filename, 'longitude');
    lon_series_2 = lon_series_2 - (lon_series_2 > 180) * 360;
    dphi_2 = abs(lat_series_2(2) - lat_series_2(1));
    dlambda_2 = abs(lon_series_2(2) - lon_series_2(1));
    [lat_indices_2, lon_indices_2] = latlonindices(lat_series_2, lon_series_2, ...
                                     latmin, latmax, lonmin, lonmax);
    lat_2 = lat_series_2(lat_indices_2);
    lon_2 = lon_series_2(lon_indices_2);

    [X, Y] = meshgrid(lon, lat);
    lon_2 = lon_2 - 360.0 * (lon_2 > 180.0);
    [X_2, Y_2] = meshgrid(lon_2, lat_2);

    % selecting precipitation events, the idea is scanning half hourly
    % precipitation data, for a certain period defined as file_period; then
    % calculate the accumulated precipitation during this period. if one
    % certain event starts to exceed the threshold at time t, define its
    % starting time as t - file_period/2.
    
    for f = 1 : length(precip_files) - file_period + 1
        if f == 1
            for ff = 1 : file_period
                filename = precip_files(f - 1 + ff).name;
                precip_temp = h5read(filename, '/Grid/HQprecipitation');
                FillValue = h5readatt(filename, '/Grid/HQprecipitation', '_FillValue');
                precip_temp = precip_temp(lat_indices, lon_indices);
                precip_temp = precip_temp .* (precip_temp ~= FillValue);
                precip_temp = interp2(X, Y, precip_temp, X_2, Y_2, 'linear');    
                precip_temp = precip_temp .* (precip_temp > 0);
                if ff == 1
                    precip = precip_temp;
                else
                    precip = precip + precip_temp;
                end
            end
        else
            filename_end = precip_files(f - 1 + file_period).name;
            precip_end = h5read(filename_end, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename_end, '/Grid/HQprecipitation', '_FillValue');
            precip_end = precip_end(lat_indices, lon_indices);
            precip_end = precip_end .* (precip_end ~= FillValue);
            precip_end = interp2(X, Y, precip_end, X_2, Y_2, 'linear', 0.0);
            precip_end = precip_end .* (precip_end > 0);
            filename_1 = precip_files(f - 1).name;
            precip_1 = h5read(filename_1, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename_1, '/Grid/HQprecipitation', '_FillValue');
            precip_1 = precip_1(lat_indices, lon_indices);
            precip_1 = precip_1 .* (precip_1 ~= FillValue);
            precip_1 = interp2(X, Y, precip_1, X_2, Y_2, 'linear', 0.0);
            precip_1 = precip_1 .* (precip_1 > 0);
            precip = precip - precip_1 + precip_end;
        end

        [lats_, lons_] = find(precip > threshold_precip);
        disp(strcat('f = ', num2str(f), ', ', num2str(length(lats_)), ' extreme point(s) found.'));
        if ~isempty(precip_events)
            for i = 1 : length(precip_events)
                precip_events(i).added = 0;
            end
        end
        while ~isempty(lats_)
            if isempty(precip_events)
                precip_events = event;
                precip_events.region_lon = lons_(1);
                precip_events.region_lat = lats_(1);
                precip_events.mean_lon = lons_(1);
                precip_events.mean_lat = lats_(1);
                precip_events.start_time = f + file_period / 2;
                precip_events.added = 1;
                disp('one precipitation event created');
                disp(strcat('[', num2str(lons_(1)), ', ', num2str(lats_(1)), ']', ...
                     ' is added to event 1'));
            else
                aaa = ([precip_events.mean_lon] - lons_(1)) * dlambda_2 / 180.0 * 3.14159 * R_earth; % in m
                bbb = ([precip_events.mean_lat] - lats_(1)) * dphi_2 / 180.0 * 3.14159 * R_earth; % in m
                mean_distance = sqrt(aaa.^2 + bbb.^2); % calculate mean distance in every event
                % find the closest event's index
                temp_index = find(mean_distance == min(mean_distance));
                % if multiple points are of the same distance, pick the first one
                temp_index = temp_index(1); 
                aaa = ([precip_events(temp_index).region_lon] - lons_(1)) * dlambda_2 / 180.0 * 3.14159 * R_earth;
                bbb = ([precip_events(temp_index).region_lat] - lats_(1)) * dphi_2 / 180.0 * 3.14159 * R_earth;
                temp_distance = sqrt(aaa.^2 + bbb.^2);
                % if the smallest distance is larger than the large-scale threshold distance, create a new event
                if min(temp_distance) > LS_distance
                    precip_events(end + 1).region_lon = lons_(1);
                    % the 'end + 1' adds a new element in the 'precip_events' array
                    precip_events(end).region_lat = lats_(1);
                    precip_events(end).mean_lon = lons_(1);
                    precip_events(end).mean_lat = lats_(1);
                    precip_events(end).start_time = f + file_period / 2;
                    precip_events(end).added = 1;
                    disp('one precipitation event created');
                    disp(strcat('[', num2str(lons_(1)), ', ', num2str(lats_(1)), ']', ...
                                ' is added to event  ', num2str(temp_index)));
                % else, if this event already contains this point
                elseif ~isempty(find((precip_events(temp_index).region_lon == lons_(1)) .* ...
                        (precip_events(temp_index).region_lat == lats_(1)), 1))
                    disp(strcat('point ', '[', num2str(lons_(1)), ', ', num2str(lats_(1)), ']', ...
                         ' already exists, ignored'));
                    precip_events(end).added = 1;
                    % skip points that are already included in the region
                % else, add this point to the closest existing event
                else
                    precip_events(temp_index).region_lon(end + 1) = lons_(1);
                    precip_events(temp_index).region_lat(end + 1) = lats_(1);
                    precip_events(temp_index).mean_lon = mean([precip_events(temp_index).region_lon]);
                    precip_events(temp_index).mean_lat = mean([precip_events(temp_index).region_lat]);
                    precip_events(end).added = 1;
                    disp(strcat('[', num2str(lons_(1)), ', ', num2str(lats_(1)), ']', ...
                         ' is added to event ', num2str(temp_index)));
                end
            end
            % after assignment, clear the array
            lats_(1) = [];
            lons_(1) = [];
        end
        if ~isempty(precip_events)
            i = 1;
            while i <= length(precip_events)
                % if this event is not added points at this time step, or comes to the end of files,
                % then count it as ended
                if (~precip_events(i).added) || f == length(precip_files) - file_period + 1
                    precip_events(i).end_time = f + file_period / 2;
                    if ~isempty(precip_events_output)
                        precip_events_output(end + 1) = precip_events(i);
                    else
                        precip_events_output = precip_events(i);
                    end
                    disp(strcat('precipitation event ', num2str(length(precip_events_output)), ' ended'));
                    precip_events(i) = [];
                    i = i - 1;
                end
                i = i + 1;
            end
        end
        %clear('precip_events');
    end

    % second, define the large-scale area of events based on a certain distance. 
    
    EXPANSION = false;
    if EXPANSION
        for k = 1 : length(precip_events_output)
            LS_index_x = ceil(LS_distance / R_earth * 180.0 / 3.14159 / dlambda_2);
            LS_index_y = ceil(LS_distance / R_earth * 180.0 / 3.14159 / dphi_2);
            up = min(max(precip_events_output(k).region_lat) + LS_index_y, length(lat_2)); 
            down = max(min(precip_events_output(k).region_lat) - LS_index_y, 1);
            left = max(min(precip_events_output(k).region_lon) - LS_index_x, 1);
            right = min(max(precip_events_output(k).region_lon) + LS_index_x, length(lon_2));
            temp_lons_ = [];
            temp_lats_ = [];
            for i = left : right
                for j = down : up
                    aaa = ([precip_events_output(k).region_lon] - i) * dlambda_2 / 180.0 * 3.14159 * R_earth;
                    bbb = ([precip_events_output(k).region_lat] - j) * dphi_2 / 180.0 * 3.14159 * R_earth;
                    temp_distance = sqrt(aaa.^2 + bbb.^2);
                    if min(temp_distance) < LS_distance
                       temp_lons_(end + 1) = i;
                       temp_lats_(end + 1) = j;
                    end
                end
            end
            precip_events_output(k).region_lon = [precip_events_output(k).region_lon, temp_lons_];
            precip_events_output(k).region_lat = [precip_events_output(k).region_lat, temp_lats_];
        end
    end
    

    % third, merge different events that has overlayed start time and region
    
    i = 1;
    while i <= length(precip_events_output) - 1
        j = i + 1;
        while j <= length(precip_events_output)
            % if event i's end time is not too far from event j's start time
            if abs(precip_events_output(i).end_time - precip_events_output(j).start_time) < ceil(file_period / 4)
                % if event j has points that is in event i
                if any(ismember(precip_events_output(i).region_lon, precip_events_output(j).region_lon) .* ...
                        ismember(precip_events_output(i).region_lat, precip_events_output(j).region_lat) == 1)
                    for k = 1 : length(precip_events_output(j).region_lon)
                        % find points that is not in event i
                        if isempty(find((precip_events_output(i).region_lon == precip_events_output(j).region_lon(k)) .* ...
                        (precip_events_output(i).region_lat == precip_events_output(j).region_lat(k)), 1))
                            precip_events_output(i).region_lon(end + 1) = precip_events_output(j).region_lon(k);
                            precip_events_output(i).region_lat(end + 1) = precip_events_output(j).region_lat(k);
                            precip_events_output(i).end_time = precip_events_output(j).end_time;
                        end
                    end
                    precip_events_output(i).mean_lon = mean([precip_events_output(i).region_lon]);
                    precip_events_output(i).mean_lat = mean([precip_events_output(i).region_lat]);
                    precip_events_output(i).end_time = precip_events_output(j).end_time;
                    precip_events_output(j) = [];
                    j = j - 1;
                end
            end
            j = j + 1;
        end
        i = i + 1;
    end
    % fourth, scan back and forth to cover the time span of this event
    
    %{
    for i = 1 : length(precip_events_output)
        f = precip_events_output(i).start_time
        while f > 1
            f = f - 1;
            filename = precip_files(f).name;
            precip_temp = h5read(filename, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename, '/Grid/HQprecipitation', '_FillValue');
            precip_temp = precip_temp(lat_indices, lon_indices);
            precip_temp = precip_temp .* (precip_temp ~= FillValue);
            precip_temp = interp2(X, Y, precip_temp, X_2, Y_2, 'linear');
            [temp_lat_, temp_lon_] = find(precip_temp(precip_events_output(i).region_lat, ...
                                     precip_events_output(i).region_lon) > threshold_precip_2);
            if length(temp_lat_) >= ratio_LS * length(precip_events_output(i).region_lat)
                precip_events_output(i).end_time = f;
            else 
                break;
            end
        end
        f = precip_events_output(i).end_time
        while f < length(precip_files)
            f = f + 1;
            filename = precip_files(f).name;
            precip_temp = h5read(filename, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename, '/Grid/HQprecipitation', '_FillValue');
            precip_temp = precip_temp(lat_indices, lon_indices);
            precip_temp = precip_temp .* (precip_temp ~= FillValue);
            precip_temp = interp2(X, Y, precip_temp, X_2, Y_2, 'linear');
            [temp_lat_, temp_lon_] = find(precip_temp(precip_events_output(i).region_lat, ...
                                     precip_events_output(i).region_lon) > threshold_precip_2);
            if length(temp_lat_) >= ratio_LS * length(precip_events_output(i).region_lat)
                precip_events_output(i).end_time = f + file_period / 2 ;
            else
                break;
            end
        end
    end
    %}
   
    % analyze the precipitation events
    

    event_period = zeros(length(precip_events_output), 1);
    for i = 1 : length(precip_events_output)
        event_period(i) = precip_events_output(i).end_time - precip_events_output(i).start_time;
    end
    event_index = find(event_period == max(event_period));
        
    PLOT = true;

    for k = 1 : length(precip_events_output)
        precip_event = precip_events_output(k);
        precip_event.accu_precip = zeros(length(precip_event.region_lon), 1);
        for f = precip_event.start_time - file_period / 2 : precip_event.end_time - 1 + file_period / 2
        %for f = precip_event.start_time: precip_event.end_time
            filename = precip_files(f).name;
            precip_temp = h5read(filename, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename, '/Grid/HQprecipitation', '_FillValue');
            precip_temp = precip_temp(lat_indices, lon_indices);
            precip_temp = precip_temp .* (precip_temp ~= FillValue);
            precip_temp = interp2(X, Y, precip_temp, X_2, Y_2, 'linear', 0.0);
            precip_temp = precip_temp .* (precip_temp > 0);
            if f == precip_event.start_time - file_period / 2
                precip = precip_temp / 2; % this is because the data is mm/h, to convert from mm/h to 
                                          % accumulated precipitation we have to have it divided by 2.
            else
                precip = precip + precip_temp / 2;
            end
            for l = 1 : length(precip_event.region_lon)
                i = precip_event.region_lon(l);
                j = precip_event.region_lat(l);
                precip_event.accu_precip(l) = precip_event.accu_precip(l) + precip_temp(j, i);
            end
        end 
        precip_events_output(k).accu_precip = ...
                    precip(sub2ind(size(precip), precip_event.region_lat, precip_event.region_lon)');
        
        % start plotting
        if PLOT == true
            disp(strcat('plotting figure ', num2str(k)));
            [Y_graph, X_graph] = meshgrid(lat_2, lon_2);
            figure('Visible','off');
            hold on;
            axesm ('eqdcylin', 'Frame', 'on', 'Grid', 'on');
            worldmap([latmin latmax], [lonmin lonmax]);
            setm(gca,'mapprojection','miller'); % set map projection type
            % note, contourfm has to be used this way, compared with standard contourf
            [c, h] = contourfm(double(lat_2), double(lon_2), double(precip), ...
                                20, 'LineStyle', 'none'); % draw contour on a map  
            cbar = colorbar;
            cbar.Label.String = 'mm';
            colormap('bone');
            colormap(flipud(colormap));
            load coastlines;
            geoshow(coastlat,coastlon, 'Color', 'black');
            plotm(double(lat_2(precip_event.region_lat)), ...
                  double(lon_2(precip_event.region_lon)), 'o', 'Color', 'red');
            caxis([0., 150.]);
            title(['Event-wise accumulated precipitation: ', ...
                    num2str(precip_event.end_time / 2 - precip_event.start_time / 2), ...
                    ' hours']);
            saveas(h, strcat('event_plot/event_', num2str(k), '_accu_precip.png'));
            hold off;			
        end
    end

    % calculate event precipitation time series

    precip_event = precip_events_output(event_index);
    f_start = precip_event.start_time - file_period / 2;
    f_end = precip_event.end_time - file_period / 2;
    f_start = 97 - file_period / 2;
    f_end = 336 - file_period / 2;
    %event_precip_series = zeros(precip_event.end_time - precip_event.start_time + 1, 1);
    event_precip_series = zeros(f_end - f_start + 1, 1);
    %for f = precip_event.start_time - file_period / 2 : precip_event.end_time - file_period / 2
    for f = f_start : f_end
        if f == f_start
            for ff = f : f + file_period - 1
                filename = precip_files(ff).name;
                precip_temp = h5read(filename, '/Grid/HQprecipitation');
                FillValue = h5readatt(filename, '/Grid/HQprecipitation', '_FillValue');
                precip_temp = precip_temp(lat_indices, lon_indices);
                precip_temp = precip_temp .* (precip_temp ~= FillValue);
                precip_temp = interp2(X, Y, precip_temp, X_2, Y_2, 'spline', 0.0);
                precip_temp = precip_temp .* (precip_temp > 0);
                if ff == f
                    precip = precip_temp / 2;
                else
                    precip = precip + precip_temp / 2; 
                end 
            end
        else
            filename_end = precip_files(f - 1 + file_period).name;
            precip_end = h5read(filename_end, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename_end, '/Grid/HQprecipitation', '_FillValue');
            precip_end = precip_end(lat_indices, lon_indices);
            precip_end = precip_end .* (precip_end ~= FillValue);
            precip_end = interp2(X, Y, precip_end, X_2, Y_2, 'spline', 0.0);
            precip_end = precip_end .* (precip_end > 0);
            filename_1 = precip_files(f - 1).name;
            precip_1 = h5read(filename_1, '/Grid/HQprecipitation');
            FillValue = h5readatt(filename_1, '/Grid/HQprecipitation', '_FillValue');
            precip_1 = precip_1(lat_indices, lon_indices);
            precip_1 = precip_1 .* (precip_1 ~= FillValue);
            precip_1 = interp2(X, Y, precip_1, X_2, Y_2, 'spline', 0.0);
            precip_1 = precip_1 .* (precip_1 > 0);
            precip = precip - precip_1 / 2 + precip_end / 2;
        end
        for l = 1 : length(precip_event.region_lon)
            i = precip_event.region_lon(l);
            j = precip_event.region_lat(l);
            event_precip_series(f - f_start + 1) = ...
            event_precip_series(f - f_start + 1) + precip(j, i); 
        end
    end
    event_precip_series = event_precip_series / length(precip_event.region_lon);


    fig = figure;
    xaxis = 1 : f_end - f_start + 1;
    hold on;
    plot(xaxis, event_precip_series, 'LineWidth',3);
    title('precipitation strength');
    %xlabel('time');
    ylabel(['normalized precipitation (mm/', num2str(file_period / 2),'hr)']);
    xticks = [0 48 96 144 192 240];
    xticklabels = {'\begin{tabular}{c}Aug 10 \\ 00:00\end{tabular}', ...
                   '\begin{tabular}{c}Aug 11 \\ 00:00\end{tabular}', ...
                   '\begin{tabular}{c}Aug 12 \\ 00:00\end{tabular}', ...
                   '\begin{tabular}{c}Aug 13 \\ 00:00\end{tabular}', ...
                   '\begin{tabular}{c}Aug 14 \\ 00:00\end{tabular}', ...
                   '\begin{tabular}{c}Aug 15 \\ 00:00\end{tabular}'};
    set(gca,'xtick', xticks, 'XTickLabel', xticklabels, 'TickLabelInterpreter', 'latex');
    set(gca,'linewidth',3)
    plot([precip_event.start_time - f_start - file_period / 2, ...
          precip_event.start_time - f_start - file_period / 2], [0, 100], '--red', 'LineWidth',3);
    plot([precip_event.end_time - f_start - file_period / 2, ...
          precip_event.end_time - f_start - file_period / 2], [0, 100], '--red', 'LineWidth',3);
    hold off;
    saveas(fig, 'time_series.png');
    
    % the event starts from Aug 11 22:30 to Aug 13 10:30  



    % read in data

    time_temp = ncread(nc_filename, 'time');
    level_temp = double(ncread(nc_filename, 'level'));
    T_temp = ncread(nc_filename, 't');
    omega_temp = ncread(nc_filename, 'w');
    T_temp = T_temp(lon_indices_2, lat_indices_2(end) : -1 : lat_indices_2(1), :, :);
    omega_temp = omega_temp(lon_indices_2, lat_indices_2(end) : -1 : lat_indices_2(1), :, :);

    omega_QG_temp = ncread('omega_QG_0.75.nc', 'omega_QG');
    omega_QG_temp = omega_QG_temp(lon_indices_2, lat_indices_2(end) : -1 : lat_indices_2(1), :, :);

    %event_index = 8;
    %precip_event = precip_events_output(event_index);

    for e = 1 : length(precip_events_output)
        precip_event = precip_events_output(e);
    
        % some selection parameters

        t_start = ceil((precip_file_start + precip_event.start_time - 1) * 0.5 / 6.0);
        t_end = floor((precip_file_start + precip_event.end_time - 1) * 0.5 / 6.0);
        if t_end < t_start
            t_end = t_start;
        end
        p_start = 1; % smaller p
        p_end = 37; % larger p

        % reverse latitide and level coordinates to satisfy right-handed xyz system
    
        lats_ = precip_event.region_lat;
        lons_ = precip_event.region_lon;
        level = level_temp(p_end : -1 : p_start) * 100; % hPa -> Pa
        phi = lat_2 / 180.0 * 3.1415926;
        lambda = lon_2 / 180.0 * 3.1415926;
        %z_ = z_temp(lons_, lats_, p_end : -1 : p_start, t_start:t_end);
        [omega_, T_, omega_QG_] = deal(zeros(length(lons_), p_end - p_start + 1, t_end - t_start + 1));
        for i = 1 : length(lons_)
            T_(i, :, :) = T_temp(lons_(i), lats_(i), p_end : -1 : p_start, t_start:t_end);
            omega_(i, :, :) = omega_temp(lons_(i), lats_(i), p_end : -1 : p_start, t_start:t_end);
            omega_QG_(i, :, :) = omega_QG_temp(lons_(i), lats_(i), p_end : -1 : p_start, t_start:t_end);
        end
        %u_ = u_temp(lons_, lats_, p_end : -1 : p_start, t_start:t_end);
        %v_ = v_temp(lons_, lats_, p_end : -1 : p_start, t_start:t_end);
        time = time_temp(t_start:t_end);

        % start calculation
    
        [Phi, Lambda] = meshgrid(phi, lambda);
    
        precip_scaling = zeros(length(lons_), t_end - t_start + 1);
        precip_scaling_QG = zeros(length(lons_), t_end - t_start + 1);
        for i = 1 : length(lons_)
            for t = 1 : length(time)
                p_t_index = 23; % set tropopause to 200hPa
                %p_t_index = 10;
                p_s_index = 3; % set surface to surface pressure
                p_index = p_s_index : p_t_index;
                p_star_i = p_star(T_(i, p_index, t));
                T_i = T_(i, p_index, t);
                r_star = gamma * p_star_i ./ level(p_index)'; % mass mixing ratio kg/kg
                A = - gamma * p_star_i ./ (gamma * p_star_i + level(p_index)').^2;
                %B1 = (R / cp ./ level(p_index)' + gamma * Lv * p_star_i ./ (cp * T_i .* level(p_index)'.^2)) ./ ...
                %     (1 ./ T_i - Lv * r_star ./ (cp * T_i.^2) + ...
                %      Lv^2 * gamma * p_star_i ./ (cp * Rv * T_i.^3 .* level(p_index)'));
                B1 = T_i * R ./ (cp * level(p_index)') .* (1 + Lv * r_star ./ (R * T_i)) ./...
                        (1 + r_star * Lv ./ (cp * T_i) .* (Lv ./ (Rv * T_i) - 1.0));
                B2 = gamma * Lv * level(p_index)' .* p_star_i ./ ...
                     ((gamma * p_star_i + level(p_index)').^2 * Rv .* T_i.^2);
                temp = omega_(i, p_s_index : p_t_index, t) .* (A + B1 .* B2);
                temp_QG = omega_QG_(i, p_s_index : p_t_index, t) .* (A + B1 .* B2);
                precip_scaling(i, t) = 0;
                precip_scaling_QG(i, t) = 0;
                for h = 1 : p_t_index - p_s_index + 1
                    if h == 1
                        delta_p = (level(h) - level(h + 1));
                    elseif h == length(level)
                        delta_p = (level(h - 1) - level(h));
                    else
                        delta_p = level(h - 1) - level(h + 1) / 2.;
                    end
                    precip_scaling(i, t) = precip_scaling(i, t) - delta_p * temp(h) / g;
                    precip_scaling_QG(i, t) = precip_scaling_QG(i, t) - delta_p * temp_QG(h) / g;
                end
            end
        end
        % 12-hourly precipitation scaling
        precip_events_output(e).precip_scaling = mean(precip_scaling, 2) * 3600 * 12;
        precip_events_output(e).precip_scaling_QG = mean(precip_scaling_QG, 2) * 3600 * 12;
        % real 12-hourly precipitation, including start - 6hr to end + 6hr
        precip_events_output(e).accu_precip / (precip_event.end_time - precip_event.start_time + 1 + file_period) * 12;
    end
    clear('lat_temp', 'z_temp', 'T_temp', 'omega_temp', 'u_temp', 'v_temp', 'time_temp');

    fig = figure;
    hold on;
    for e = 1 : length(precip_events_output)
        scatter(mean(precip_events_output(e).accu_precip), mean(precip_events_output(e).precip_scaling));
    end
    title('accumulated precipitation vs. scaling');
    xlabel([num2str(file_period / 2), 'hr accumulated precipitation (mm)']);
    ylabel([num2str(file_period / 2), 'hr precipitation scaling (mm)']);
    saveas(fig, 'scaling.png');
    hold off;

	fig = figure;
    hold on;
    for e = 1 : length(precip_events_output)
        scatter(mean(precip_events_output(e).accu_precip), mean(precip_events_output(e).precip_scaling_QG));
    end
    title('accumulated precipitation vs. QG scaling');
    xlabel([num2str(file_period / 2), 'hr accumulated precipitation (mm)']);
    ylabel([num2str(file_period / 2), 'hr QG precipitation scaling (mm)']);
    saveas(fig, 'scaling_QG.png');
    hold off;

%end

function p_s = p_star(T)
        
    Rv = 462.48;
    Lv = 24.0e5;
    p_s0 = 101325.0;
    T0 = 373.0;
    p_s = p_s0 * exp(Lv / Rv * (1 / T0 - 1 ./ T));

end



