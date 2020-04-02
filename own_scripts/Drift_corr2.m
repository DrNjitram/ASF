%load relevant data file
load('res_fov1.mat')
% create useful matrices
allDeltas = zeros(1,4);
res_yes = [res(:,1),res(:,2),res(:,6),res(:,8)];

% for every bead (loop 1), calculate movement between every frame (loop 2)
for i = 1:1873 %1873 is the last bead ID
    vector = res_yes(res_yes(:,4)==i);
    FOV = length(vector);
    deltaXY = zeros(FOV-1,3);
    for j = 1:FOV-1
        deltaXY(j,1:2) = res_yes(j+1,1:2)- res_yes(j,1:2);
        deltaXY(j,3) = j; % adds the timestep (frame j+1 minus framej)
        testDelta{j} = deltaXY(j,:);
    end
    beadID = ones(FOV-1,1)*i;
    deltaXY = [deltaXY, beadID]; % adds the beadID 
    allDeltas = [allDeltas; deltaXY];
    clear deltaXY
    disp(i/1873)
end
disp('done')

