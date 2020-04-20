clear variables
clear figures
%add paths
addpath('own_scripts');  

% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = [pwd '\data\'];
addpath([pwd '\Code_2D_feature_finding\']);
featsize = 3; 
masscut = 3000;
create_intermediate_graphs = true;

frame = 1;
fovn = FOV;

% These variables can be used to check the quality using %[M2,MT] =
% mpretrack_init(basepath, featsize, barint, barrg, barcc, IdivRg, fovn,
% frame);
barint = 4000;
barrg = 4;
barcc = 0.5;
IdivRg = barint/barrg;

% Define filename
fname = [basepath 'colloids.tif'];

% Get tiff information
disp("Reading Tiff");
info = imfinfo(fname);
num_images = numel(info);
fprintf('Found %d frames\n\n',num_images);
% Uncomment below line to only examine the first image
%num_images = 1;

% Read tiff stack into memory

images = readtiffstack(fname);

% Create cell to contain features
r = {num_images};

% For each image, get the features
fprintf('Finding features:\nProcessing frame ');
linelength = 0;
for i = 1:num_images
    fprintf(repmat('\b',1,linelength));
    linelength = fprintf('%d of %d', i, num_images);
    
    r{i} = feature2D(images{i}, fovn, featsize, masscut, barrI);
end
fprintf('\n\n');

%% Construct MT for further processing
MT = [];

for i = 1:num_images
    MT = [MT; r{i}, ones(length(r{i}(:,1)), 1)*i, ones(length(r{i}(:,1)), 1)*frametime*(i)];
end

% Save created matrix
if ~exist('data/Feature_finding', 'dir'); mkdir('data/Feature_finding'); end
addpath('data/Feature_finding');
save([basepath 'Feature_finding/MT_' num2str(FOV) '_Feat_Size_' num2str(featsize) '.mat'], 'MT');

%% Tracking
fancytrack(basepath, FOV, featsize, 2, 70, 1);

%% Calculate standard deviation of remainders mod 1
remainder_x = [];
remainder_y = [];

for i = 1:num_images
    remainder_x = [remainder_x, mod(r{i}(:,1), 1)'];
    remainder_y = [remainder_y, mod(r{i}(:,2), 1)'];
end

%% Imaging
if create_intermediate_graphs
    % Frame to display
    no = 1;
    figure
    imshow(images{i}, [0, max(max(images{i}))]);
    hold on
    plot(r{no}(:,1), r{no}(:,2), '.');

    % Display histograms of remainders
    figure
    histogram(remainder_x, 10)
    figure
    histogram(remainder_y, 10)
end



