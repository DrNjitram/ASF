clear 

% Define filename
fname = 'data\colloids.tif';

% Get tiff information
info = imfinfo(fname);
num_images = numel(info);
% Uncomment below line to only examine the first image
%num_images = 1;

% Read tiff stack into memory
images = readtiffstack(fname);

% Create cell to contain features
r = {num_images};

% For each image, get the features
for i = 1:num_images
    r{i} = feature2D(images(:, :, i),1,3,3000);
end

%% Calculate standard deviation of remainders mod 1
remainder_x = [];
remainder_y = [];

for i = 1:num_images
    remainder_x = [remainder_x, std(mod(r{i}(:,1), 1)')];
    remainder_y = [remainder_y, std(mod(r{i}(:,2), 1)')];
end

%% Imaging

% Frame to display
no = 1;
figure
imshow(images(:, :, no), [0, max(max(images(:, :, no)))]);
hold on
plot(r{no}(:,1), r{no}(:,2), '.');

% Display histograms of remainders
figure
histogram(remainder_x, 10)
figure
histogram(remainder_y, 10)