clear variables

% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = 'D:\Uni\ASM\data\';
frametimes = 4.2; %s
numframes = 79;

time = [1:79]*frametimes;
save([basepath 'fov' num2str(FOV) '_times.mat'], 'time');

featsize = 3; 
Imin = 3000;
frame = 1;
fovn = FOV;
barint = 4000;
barrg = 4;
barcc = 0.5;
IdivRg = barint/barrg;

%[M2,MT] = mpretrack_init(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame);


mpretrack( basepath, fovn, featsize, barint, barrg, barcc, IdivRg, numframes);

fancytrack(basepath, FOV, featsize);