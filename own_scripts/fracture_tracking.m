clear variables
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
    
    r{i} = feature2D(images{i}, fovn, featsize, masscut);
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
    remainder_x = [remainder_x, std(mod(r{i}(:,1), 1)')];
    remainder_y = [remainder_y, std(mod(r{i}(:,2), 1)')];
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
%% Driftcorrection

% INPUT
% 
% matrix called res, originating from 'fancytrack.m' ADD MORE DETAILS
% CONCERNING MATRIX FORMAT

% OUTPUT
% 
% 


load('res_fov1.mat') %load relevant data file
res_yes = [res(:,1),res(:,2),res(:,6),res(:,8)]; %remove irrelevant values

%% drift correction based on average positions
% INPUT: [res_yes]
% OUTPUT: []

sortFOV = sortrows(res_yes, 3); %sorts the matrix to allow easy selection of frame i
corMatrix = sortFOV(:,1:2); %create array to hold subtracted averages
corMatrix(:,:) = 1;
corMatrix(1,:) = 0;
maxFOV = max(sortFOV(:,3)); %gives final frame for defining for loop below
for i = 1:maxFOV
    FOV = find(sortFOV(:,3) == i); %gives vector containing every row belonging to frame i
    maxRow = max(FOV);
    minRow = min(FOV);
    X = sortFOV(minRow:maxRow,1); %extracts x-values belonging to frame i
    Y = sortFOV(minRow:maxRow,2); %same as above for y-values
    if i ~= 1 %nothing to be subtracted from first frame
        meanX(i) = mean(X); %calculates average x-position
        meanY(i) = mean(Y); %calculates average Y-position
        corMeanX = meanX(i)-meanX(i-1); %subtracts average position of the frame before from the average position of a frame.
        corMeanY = meanY(i)-meanY(i-1); %same as above but for Y
    end
    corMatrix(minRow:maxRow,1) = corMatrix(minRow:maxRow,1)*corMeanX; %puts the x-values we use for correction at the correct place
    corMatrix(minRow:maxRow,2) = corMatrix(minRow:maxRow,2)*corMeanY; %same as above for y-values
end
corSortFOV = [sortFOV(:,1:2)- corMatrix(:,1:2),sortFOV(:,3:4)]; %subtracts the correction value from the original positions to give the drift corrected positions
corRes_yes = sortrows(corSortFOV,4); % sorted by beadID again

%% Next step: Show the movement plots of particles (Tom's code)
% Plots 50 particles as a function of frame for both pre-driftcorrection
% (figure 1) and post-driftcorrection (figure 2)

figure(1) %pre-correction
for m=1:50;
    plot(res((1+(79*(m-1))):(79*m),1),res((1+(79*(m-1))):(79*m),2)); hold on;
end
figure(2) %post-correction
for m=1:50;
    plot(corRes_yes((1+(79*(m-1))):(79*m),1),corRes_yes((1+(79*(m-1))):(79*m),2)); hold on;
end

%% drift correction failures (remove if agreed upon)
% Selecting beadIDs The first loop goes through all beadIDs and creates a
% new array  per bead containing only the positions for the relevant bead
% and adds data calculated in the second loop to allData.
% for i = 1:max(res_yes(:,4))
%     rows = find(res_yes(:,4)==i);
%     FOV = length(rows);
%     lastRow = max(rows);
%     firstRow = min(rows);
%     XY = res_yes(firstRow:lastRow,1:2);
%     deltaXY = zeros(FOV-1,3);
%     % Selecting all frames for current beadID The second loop calculates
%     % the displacement between every subsequent frame.
%     for j = 1:length(XY)-1 
%         deltaXY(j,1:2) = XY(j+1,1:2) - XY(j,1:2);
%         deltaXY(j,3) = j; % adds dFOV
%     end
%     beadID = ones(FOV-1,1)*i;
%     deltaXY = [deltaXY, beadID]; % adds the beadID 
%     allDeltas = [allDeltas; deltaXY];
%     clear deltaXY
%     disp(i/max(res_yes(:,4)))
% end
% disp('done')
%%
% The next step is to decide on a(n average) value per dFOV to subtract
% from the the original positions INPUT: [allDeltas] OUTPUT: []
% 
% sortdFOV = sortrows(allDeltas, 3);
% corValues = sortdFOV;
% maxdFOV = max(sortdFOV(:,3));
% for i = 1:maxdFOV
%     dFOV = find(sortdFOV(:,3) == i);
%     maxRow = max(dFOV);
%     minRow = min(dFOV);
%     averagedX(i) = mean(sortdFOV(minRow:maxRow,1));
%     averagedY(i) = mean(sortdFOV(minRow:maxRow,2));
%     corValues(minRow:maxRow,1) = sortdFOV(minRow:maxRow,1) - averagedX(i);
%     corValues(minRow:maxRow,2) = sortdFOV(minRow:maxRow,2) - averagedY(i);
% end
% corValues(:,3:4) = sortdFOV(:,3:4);
% corAllDeltas = sortrows(corValues,4);