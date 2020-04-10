clear variables
close all

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
corMatrix = zeros(length(sortFOV), 2); %create array to hold subtracted averages

maxFOV = max(sortFOV(:,3)); %gives final frame for defining for loop below

meanX = zeros(1, maxFOV);
meanY = zeros(1, maxFOV);
corMeanX = zeros(1, maxFOV);
corMeanY = zeros(1, maxFOV);

for i = 1:maxFOV
    FOV = find(sortFOV(:,3) == i); %gives vector containing every row belonging to frame i
    maxRow = max(FOV);
    minRow = min(FOV);
    X = sortFOV(minRow:maxRow,1); %extracts x-values belonging to frame i
    Y = sortFOV(minRow:maxRow,2); %same as above for y-values

    meanX(i) = mean(X); %calculates average x-position
    meanY(i) = mean(Y); %calculates average Y-position
    
    if i == 1
        corMeanX(1) = 0;
        corMeanY(1) = 0;
    else
        corMeanX(i) = meanX(i)-meanX(i-1); %subtracts average position of the frame before from the average position of a frame.
        corMeanY(i) = meanY(i)-meanY(i-1); %same as above but for Y
    end
    
    corMatrix(minRow:maxRow,1) = sum(corMeanX(1:i)); %puts the x-values we use for correction at the correct place
    corMatrix(minRow:maxRow,2) = sum(corMeanY(1:i)); %same as above for y-values
end

corSortFOV = [sortFOV(:,1:2)- corMatrix(:,1:2),sortFOV(:,3:4)]; %subtracts the correction value from the original positions to give the drift corrected positions
corRes_yes = sortrows(corSortFOV,4); % sorted by beadID again

%% Next step: Show the movement plots of particles (Tom's code)
% Plots 50 particles as a function of frame for both pre-driftcorrection
% (figure 1) and post-driftcorrection (figure 2)
figure
plot([1:maxFOV], corMeanX); hold on;
plot([1:maxFOV], corMeanY); hold on;


figure %pre-correction
p = 5;% particles to plot
for m=1:p
    FOV = find(res(:,8) == m);
    plot(res(min(FOV):max(FOV), 1), res(min(FOV):max(FOV), 2)); hold on;
end
figure %post-correction
for m=1:p
    FOV = find(res(:,8) == m);
    plot(corRes_yes(min(FOV):max(FOV), 1), corRes_yes(min(FOV):max(FOV), 2)); hold on;
end