clear variables

% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = [pwd '\data\'];
addpath([pwd '\Code_2D_feature_finding\']);
frametimes = 4.2; %s
numframes = 79;

time = [0:numframes-1]*frametimes;
save([basepath 'fov' num2str(FOV) '_times.mat'], 'time');
 
masscut = 3000;
frame = 1;
fovn = FOV;
barrg = 4;
barcc = 0.5;
featsize = 3;
barrI = 1000;
Imin = 100;
IdivRg = barrI/barrg;
maxdisp = 2;
goodenough = 79;
memory = 1;

%load('locating_variables.mat')

mpretrack(basepath, fovn, featsize, barrI, barrg, barcc, IdivRg, numframes, masscut, Imin, 2);


load([basepath '\Feature_finding\MT_1_Feat_Size_',num2str(floor(featsize)),'.mat']) %concatenate strings using featsize variable
frames = find(MT(:,6)==2); %find where (which element) where the second frame starts

%%
%check visually if particle locations match the first frame
first_image = imread([basepath '\fov1\fov1_0001.tif']);
imagesc(first_image); hold on;
scatter(MT(1:(frames(1)-1),1),MT(1:(frames(1)-1),2),'*k'); 


%%
%check pixel biasing
xbias = mod(MT(:,1),1); %find using modulus, find the remainder
ybias = mod(MT(:,2),1); 

figure
histogram(xbias,10,'FaceColor','r'); hold on;
histogram(ybias,10,'FaceColor','m');
%%
% tracking 
fancytrack( basepath, fovn, featsize, maxdisp, goodenough, memory ) %all inputs required


%%
load('D:\Uni\ASM\data\Bead_tracking\res_files\res_fov1.mat')
% Plot 50 of those particles
figure
for m=1:50
    plot(res((1+(79*(m-1))):(79*m),1),res((1+(79*(m-1))):(79*m),2)); hold on;
end


%%
% dedrifting code
% explained poorly in description document; take the average position of
% all particles at t=2, and subtract the previous average position at t=1;
% this will give the displacement due to drift. Subtract this displacement
% from each particle position at t=2. 

%%
% calculate local strain for each particle