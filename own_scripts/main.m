%clear variables
close all

% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = [pwd '\data\'];
addpath([pwd '\Code_2D_feature_finding\']);
frametimes = 4.2; %s
numframes = 79;

time = [0:numframes-1]*frametimes;
save([basepath 'fov' num2str(FOV) '_times.mat'], 'time');

featsize = 3; 
masscut = 3000;
frame = 1;
fovn = FOV;
barint = 4000;
barrg = 4;
barcc = 0.5;
IdivRg = barint/barrg;

%[M2,MT] = mpretrack_init(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame);


mpretrack(basepath, fovn, featsize, barint, barrg, barcc, IdivRg, numframes, masscut);

fancytrack(basepath, FOV, featsize, 2, 70, 1);

dedrifted = dedrift(false);

final = Strain_calc(dedrifted);