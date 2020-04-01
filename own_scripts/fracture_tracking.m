clear variables
%add paths
addpath('Code_2D_feature_finding')  
addpath('own_scripts')  

% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = 'D:\Uni\ASM\';
featsize = 3; 
Imin = 3000;
create_intermediate_graphs = true;

% Define filename
fname = 'data\colloids.tif';

% Get tiff information
disp("Reading Tiff");
info = imfinfo(fname);
num_images = numel(info);
fprintf('Found %d frames\n',num_images);
% Uncomment below line to only examine the first image
%num_images = 1;

% Read tiff stack into memory

images = readtiffstack(fname);

% Create cell to contain features
r = {num_images};

% For each image, get the features
fprintf('Processing images:\nProcessing frame ');
linelength = 0;
for i = 1:num_images
    fprintf(repmat('\b',1,linelength));
    linelength = fprintf('%d of %d', i, num_images);
    
    r{i} = feature2D(images(:, :, i),1,featsize,Imin);
end
fprintf('\n');

%% Construct MT for further processing
MT = [];

for i = 1:num_images
    MT = [MT; r{i}(:,1), r{i}(:,2), ones(length(r{i}(:,1)), 1)*frametime*(i-1)];
end

% Save created matrix 
if ~exist('Feature_finding', 'dir'); mkdir('Feature_finding'); end
addpath('Feature_finding');
save([basepath 'Feature_finding/MT_' num2str(FOV) '_Feat_Size_' num2str(featsize) '.mat'], 'MT');

%% Tracking
fancytrack(basepath, FOV, featsize);

%% Calculate standard deviation of remainders mod 1
remainder_x = [];
remainder_y = [];

for i = 1:num_images
    remainder_x = [remainder_x, std(mod(r{i}(:,1), 1)')];
    remainder_y = [remainder_y, std(mod(r{i}(:,2), 1)')];
end

%% Imaging
if create_intermediate_graphs
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
end