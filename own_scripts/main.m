clear variables
close all

%% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = [pwd '\data\'];
addpath([pwd '\Code_2D_feature_finding\']);
start_frame = 10;
end_frame = 79;
numframes = end_frame - start_frame + 1;
frametimes = 4.2; %s
time = [0:numframes-1]*frametimes;

save([basepath 'fov' num2str(FOV) '_times.mat'], 'time');

featsize = 3; 
masscut = 3000;
frame = 1;
fovn = FOV;
barrI = 1000;
barrg = 4;
barcc = 0.5;
IdivRg = barrI/barrg;
maxdisp = 4; % def 1
goodenough = numframes - 2;
memory = 3; % def 2

p_size = 5;


%%
%[M2,MT] = mpretrack_init(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame);

mpretrack(basepath, start_frame, fovn, featsize, barrI, barrg, barcc, IdivRg, numframes, masscut);


fancytrack(basepath, FOV, featsize, maxdisp, goodenough, memory);

%%
dedrifted = dedrift(basepath, false);

final = Strain_calc(dedrifted, p_size, true);


%% plotting

load([basepath '\Feature_finding\MT_1_Feat_Size_',num2str(floor(featsize)),'.mat']) %concatenate strings using featsize variable
load( [basepath 'Bead_tracking/res_files/res_fov1.mat'] ); %load relevant data file

%check visually if particle locations match the first frame
frame_to_check = 60;
frames = find(MT(:,6)==frame_to_check-start_frame+1); %find where (which element) where the second frame starts
first_image = imread([basepath '\fov1\fov1_' num2str(frame_to_check,'%04i') '.tif']);
imagesc(first_image); hold on;
scatter(MT((frames),1), MT((frames),2),'*k'); 

pixel_bias(res);

n = 5; % frame number
figure
scatter(final(:,1,n), -1*final(:,2,n), 20 , final(:,4,n),'filled');
colorbar;
caxis([-0.05, 0.05]);
% load([basepath 'Feature_finding/MT_' num2str(fovn) '_Feat_Size_' num2str(featsize) '.mat'])
% load('res_fov1.mat') %load relevant data file
% 
% sortFOV = sortrows(MT, 6);
% rows = find(sortFOV(:,6) == n); %gives vector containing every row belonging to frame i
% X = sortFOV(min(rows):max(rows),1); %extracts x-values belonging to frame i
% Y = sortFOV(min(rows):max(rows),2); %same as above for y-values
% 
% figure
% plot(X, -Y, '.'); 
% 
% figure
% 
% sortFOV = sortrows(res, 6);
% 
% rows = find(sortFOV(:,6) == n); %gives vector containing every row belonging to frame i
% X = sortFOV(min(rows):max(rows),1); %extracts x-values belonging to frame i
% Y = sortFOV(min(rows):max(rows),2); %same as above for y-values
% plot(X, -Y, '.'); 






%% Plot process zone for frames of interest
final = final;  %the xx xy yy yx strain matrix
t = 7;          %the frame you want to see
strain_type = 3;%choose 3/4/5/6: 3 = xx; 4 = xy; 5 = yy 6 = yx
strain_cut = 0; %choose from below which strain value you want to exclude particles

calc_processZone(final,t,strain_type,strain_cut);

% create gif/video of all frames(if time allows)
% for i = time   