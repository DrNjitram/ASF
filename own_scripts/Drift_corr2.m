% SHOULD MAKE FUNCTION OF THIS

% INPUT
% 
% matrix called res, originating from 'fancytrack.m'
% ADD MORE DETAILS CONCERNING MATRIX FORMAT

OUTPUT

% The function outputs an array called 'allDeltas' which contains the
% displacement between frames for every beadID
% ADD MORE DETAILS CONCERNING MATRIX FORMAT
% 
% 

%load relevant data file
load('res_fov1.mat')
% create useful matrices
allDeltas = zeros(1,4);
res_yes = [res(:,1),res(:,2),res(:,6),res(:,8)];

%% Selecting beadIDs
% The first loop goes through all beadIDs and creates a new array  per bead
% containing only the positions for the relevant bead and adds data 
% calculated in the second loop to allData.
for i = 1:max(res_yes(:,4))
    rows = find(res_yes(:,4)==i);
    FOV = length(rows);
    lastRow = max(rows);
    firstRow = min(rows);
    XY = res_yes(firstRow:lastRow,1:2);
    deltaXY = zeros(FOV-1,3);
    %% Selecting all frames for current beadID
    % The second loop calculates the displacement between every subsequent
    % frame.
    for j = 1:length(XY)-1 
        deltaXY(j,1:2) = XY(j+1,1:2) - XY(j,1:2);
        deltaXY(j,3) = j; % adds the timestep (frame j+1 minus framej)
        testDelta = deltaXY(j,:);
    end
    beadID = ones(FOV-1,1)*i;
    deltaXY = [deltaXY, beadID]; % adds the beadID 
    allDeltas = [allDeltas; deltaXY];
    clear deltaXY
    disp(i/max(res_yes(:,4)))
end
disp('done')

