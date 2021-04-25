
addpath('source/')
plot_path = 'plots/';
if ~exist(plot_path)
    mkdir(plot_path);
end

downscale_ratio = 0.5;
ub = 80;
lb = 20;
r = 500 * downscale_ratio;
r_particle = 3;
numcols = 1920 * downscale_ratio;
numrows = 1090 * downscale_ratio;
inner = 1/3.5;
outer = 0.8;

v = VideoReader('experiement.m4v');
vidHeight = v.Height;
vidWidth = v.Width;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

k = 1;
interval = 0.25;
for t = 10 : interval : 180
    if mod(t, 10) == 0
        disp(['t = ', num2str(t)])
    end
    v.CurrentTime = t; % in seconds
    temp_image = rgb2gray(readFrame(v));
    s(k).cdata = imresize(temp_image, [numrows numcols]);
    k = k+1;
end

x = [1 : length(s(1).cdata(1, :))]';
y = [1 : length(s(1).cdata(:, 1))]';

x0_d = 937/1920*length(x);
y0_d = 520/1090*length(y);
x0_m = downscale_ratio*(length(x) + 1);
y0_m = downscale_ratio*(length(y) + 1);

[X, Y] = meshgrid(x, y);

% define an image for calibration
image0 = s(1).cdata;
[x0_0, y0_0] = find_center(X, Y, x0_d, y0_d, r, image0);
image0_masked = masking(x0_0, y0_0, X, Y, r, image0, inner, outer);
image0_shifted = uint8(interp2(X, Y, double(image0_masked), X + (x0_0 - x0_m), Y + (y0_0 - y0_m), 'cubic', 0));

theta0 = 96.5;
image0_final = linear_masking(downscale_ratio, theta0, X, Y, image0_shifted);
%figure
%imshow(image0_shifted);
figure
imshow(image0_final);


particles = {};

for i = 1 : length(s)
    i
    temp_image = s(i).cdata;
    
    % 1. find the center of the image
    [x0, y0] = find_center(X, Y, x0_d, y0_d, r, temp_image);
    
    % 2. mask the image to prepare for tracking
    image_masked = masking(x0, y0, X, Y, r, temp_image, inner, outer);
    
    % 3. shift the image to the center of the graph
    image_shifted = uint8(interp2(X, Y, double(image_masked), X + (x0 - x0_m), Y + (y0 - y0_m), 'cubic', 0));
    
    % 4. find the orientation of the image
    theta = principle_angle(X, Y, x0_m, y0_m, r, image_shifted, image0_final, 'sample');
    
    % 5. rotate the images
    image_rotated = imrotate(uint8(image_shifted), theta, 'crop');
    
    % 6. maskout the metal bar
    image_final = linear_masking(downscale_ratio, theta0, X, Y, image_rotated);
    s(i).image_final = image_final;
    
    % 7. find the particles
    particles{i} = find_particles(image_final, x0_m, y0_m, ub, lb, r_particle);
    
end

set(0, 'DefaultFigureVisible', 'off')
% find tracks
% 1. initiate track in the first frame
tracks = {};
working_tracks = {};
Pos = [particles{1}(:).pos];
for p = 1 : length(Pos(1, :))
    working_tracks{p} = track;
    working_tracks{p}.x_series(end + 1) = Pos(1, p);
    working_tracks{p}.y_series(end + 1) = Pos(2, p);
    working_tracks{p}.tail = [Pos(1, p); Pos(2, p)];
end
figure
hold on
contour(s(1).image_final)
scatter(Pos(1, :), Pos(2, :), 'red')
hold off
saveas(gca, [plot_path, num2str(1, '%.3d')], 'png')
% 2. append track from the following frames
set(0, 'DefaultFigureVisible', 'off')
for i = 2 : length(particles)
    
    Pos = [particles{i}(:).pos];
    
    figure
    hold on
    contour(s(i).image_final)
    scatter(Pos(1, :), Pos(2, :), 'red')
    hold off
    saveas(gca, ['plots/', num2str(i, '%.3d')], 'png')
    
    Tails = [];
    for j = 1 : length(working_tracks)
        Tails = [Tails, working_tracks{j}.tail];
    end
    added_working_tracks = {};
    for p = 1 : length(Pos(1, :))
        if ~isempty(working_tracks)
            dis = sum((repmat([Pos(1, p); Pos(2, p)], 1, length(working_tracks)) - Tails).^2, 1);
            ind = find(dis == min(dis));
            if dis(ind) < 10.0 * r_particle^2
                working_tracks{ind}.x_series(end + 1) = Pos(1, p);
                working_tracks{ind}.y_series(end + 1) = Pos(2, p);
                working_tracks{ind}.tail = [Pos(1, p); Pos(2, p)];
                Tails(:, ind) = NaN;
            else
                added_working_tracks{end + 1} = track;
                added_working_tracks{end}.x_series(end + 1) = Pos(1, p);
                added_working_tracks{end}.y_series(end + 1) = Pos(2, p);
                added_working_tracks{end}.tail = [Pos(1, p); Pos(2, p)];
            end
        else
            added_working_tracks{end + 1} = track;
            added_working_tracks{end}.x_series(end + 1) = Pos(1, p);
            added_working_tracks{end}.y_series(end + 1) = Pos(2, p);
            added_working_tracks{end}.tail = [Pos(1, p); Pos(2, p)];
        end
        
    end
    tracks = [tracks, working_tracks(~isnan(Tails(1, :)))];
    %if ~isempty(tracks)
    %    break
    %end
    working_tracks = [added_working_tracks, working_tracks(isnan(Tails(1, :)))];
end
tracks = [tracks, working_tracks];

% 3. filter the tracks
% if a track is too short, or composed of too few points, remove it
tracks_reserve = tracks;
i = 1;
while i <= length(tracks)
    temp_x = tracks{i}.x_series;
    temp_y = tracks{i}.y_series;
    dis = (temp_x' - temp_x) .* (temp_x' - temp_x) + (temp_y' - temp_y) .* (temp_y' - temp_y);
    if max(dis(:)) < 10 * r_particle^2 || length(temp_x) < 6
        tracks = {tracks{1:i-1}, tracks{i+1:end}};
    else
        i = i + 1;
    end
end

% draw the selected tracks

set(0, 'DefaultFigureVisible', 'on')
figure
for i = 1 : length(tracks)
    [theta, rho] = cart2pol(tracks{i}.x_series - x0_m, tracks{i}.y_series - y0_m);
    p = polarplot(theta, rho);
    if i == 1
        hold on
    end
    tracks{i}.theta = theta;
    tracks{i}.rho = rho;
end
hold off
saveas(gca, [plot_path, 'tracks_polar'], 'png')

figure
for i = 1 : length(tracks)
    plot(tracks{i}.x_series - x0_m, tracks{i}.y_series - y0_m)
    if i == 1
        hold on
    end
end
hold off
saveas(gca, [plot_path, 'tracks'], 'png')

% 4. transform from the tracks to speed

v_theta = [];
for i = 1 : length(tracks)
    dtheta = tracks{i}.theta(3 : end) - tracks{i}.theta(1 : end - 2);
    dtheta = mod(dtheta, 2 * pi);
    dtheta(dtheta > 6) = dtheta(dtheta > 6) - 2 * pi;
    drho = tracks{i}.rho(3 : end) - tracks{i}.rho(1 : end - 2);
    temp_v = dtheta / (2 * interval) .* tracks{i}.rho(2 : end - 1);
    %if max(abs(temp_v)) > 400
    %    break
    %end
    v_theta = [v_theta; [tracks{i}.rho(2 : end - 1)', temp_v']];
end

nbins = 50;
[N, edges, bin] = histcounts(v_theta(:, 1),nbins);
[v, v_std] = deal(zeros(nbins, 1));
for b = 1 : nbins
    v(b) = mean(v_theta(bin == b, 2));
    v_std(b) = std(v_theta(bin == b, 2));
end

save('main.mat');

figure
errorbar((edges(2:end) + edges(1:end-1))/2, v, v_std)
xlabel('radius (pixels)')
ylabel('velocity (pixels/s)')
saveas(gca, [plot_path, 'velocities'], 'png')


