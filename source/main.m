
v = VideoReader('1027_converted.m4v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first, interpolate and calibrate the image in a polar coordinate 

%while hasFrame(v)
%    video = readFrame(v);
%end
%whos video
downscale_ratio = 0.5;
ub = 80;
lb = 20;
r_particle = 3;
numcols = 1920 * downscale_ratio;
numrows = 1090 * downscale_ratio;
image0 = rgb2gray(read(v, 500));
%image1 = rgb2gray(read(v, 4000));
%image2 = rgb2gray(read(v, 4100));
image1 = rgb2gray(read(v, 30000));
image2 = rgb2gray(read(v, 31000));
image0 = imresize(image0, [numrows numcols]);
image1 = imresize(image1, [numrows numcols]);
image2 = imresize(image2, [numrows numcols]);
image0d = double(image0);
image1d = double(image1);
image2d = double(image2);
x = [1 : length(image1(1, :))]';
y = [1 : length(image1(:, 1))]';

x0_d = 937/1920*length(x);
y0_d = 520/1090*length(y);
x0_m = downscale_ratio*(length(x) + 1);
y0_m = downscale_ratio*(length(y) + 1);
theta0 = 0.0;
r = 500 * downscale_ratio;

[X, Y] = meshgrid(x, y);
%[Gmag1, Gdir1] = imgradient(image1);
%[Gmag2, Gdir2] = imgradient(image2);

% 1. find the center of the image
[x0_0, y0_0] = find_center(X, Y, x0_d, y0_d, r, image0);
[x0_1, y0_1] = find_center(X, Y, x0_d, y0_d, r, image1);
[x0_2, y0_2] = find_center(X, Y, x0_d, y0_d, r, image2);

% 2. mask the image to prepare for tracking
inner = 1/3.5;
outer = 0.8;
image0_masked = masking(x0_0, y0_0, X, Y, r, image0, inner, outer);
image1_masked = masking(x0_1, y0_1, X, Y, r, image1, inner, outer);
image2_masked = masking(x0_2, y0_2, X, Y, r, image2, inner, outer);

% 3. shift and rotate the image

image0_shifted = uint8(interp2(X, Y, double(image0_masked), X + (x0_0 - x0_m), Y + (y0_0 - y0_m), 'cubic', 0));
%contour(X - (x0_1 - x0_d), Y - (y0_1 - y0_d), image1_final)
image1_shifted = uint8(interp2(X, Y, double(image1_masked), X + (x0_1 - x0_m), Y + (y0_1 - y0_m), 'cubic', 0));
%contour(X - (x0_1 - x0_d), Y - (y0_1 - y0_d), image1_final)
image2_shifted = uint8(interp2(X, Y, double(image2_masked), X + (x0_2 - x0_m), Y + (y0_2 - y0_m), 'cubic', 0));
%contour(X - (x0_1 - x0_d), Y - (y0_1 - y0_d), image2_final)

% 4. find the orientation of the image
theta1 = principle_angle(X, Y, x0_1, y0_1, r, image1_shifted, image0_shifted, 'sample');
theta2 = principle_angle(X, Y, x0_2, y0_2, r, image2_shifted, image0_shifted, 'sample');

% 5. rotate the images
image1_final = imrotate(uint8(image1_shifted), (theta1 - theta0), 'crop');
image2_final = imrotate(uint8(image2_shifted), (theta2 - theta0), 'crop');

% compare what is before rotation and shifting and what is after these
% manipulations
figure
contour(abs(image2_final - image1_final))
figure
contour(image2 - image1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second, track particles

% 1. mask out the where the metal bar is located
principle_angle(X, Y, x0_m, y0_m, r, image1_final, image0_shifted, 'sample')
principle_angle(X, Y, x0_m, y0_m, r, image2_final, image0_shifted, 'sample')
temp1 = linear_masking(downscale_ratio, 135, X, Y, image1_final);
temp2 = linear_masking(downscale_ratio, 135, X, Y, image2_final);

%figure;
%contour(temp1<100)
figure;
imshow(temp1)
figure;
imshow(temp2)

temp1 = masking(x0_m, y0_m, X, Y, r, temp1, inner*1.01, outer*0.99);
temp2 = masking(x0_m, y0_m, X, Y, r, temp2, inner*1.01, outer*0.99);
particles1 = find_particles(temp1, x0_m, y0_m, ub, lb, r_particle);
particles2 = find_particles(temp2, x0_m, y0_m, ub, lb, r_particle);
figure
contour(temp1)
hold on
for i = 1 : length(particles)
    scatter(particles1(i).pos(1), particles1(i).pos(2))
end
hold off

[theta,rho] = cart2pol(X, Y);

temp_ind = ~ind;
Gdir1(temp_ind) = NaN;
%GG = sort([Gdir1(:), Gmag1(:)], 1);
[N, edges] = histcounts(Gdir1(:), 100);
Gx = Gmag1 .* cos(Gdir1);
Gy = Gmag1 .* sin(Gdir1);
%L1 = smooth2a(del2(image1, 0.25), 3, 3);
L2 = smooth2a(del2(image2, 0.25), 3, 3);

contourf(Gmag1 .* ind_1, 'linestyle', 'none')



