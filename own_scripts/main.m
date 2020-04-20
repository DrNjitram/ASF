clear variables
close all

%% Global variables
frametime = 4.2; %s
FOV = 1; % Unknown 
basepath = [pwd '\data\'];
addpath([pwd '\Code_2D_feature_finding\']);
start_frame = 20;
end_frame = 50;
numframes = end_frame - start_frame + 1;
time = [0:numframes-1]*frametime;

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
goodenough = numframes - 5;
memory = 4; % def 2

p_size = 20;


%%
%[M2,MT] = mpretrack_init(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame);

mpretrack(basepath, start_frame, fovn, featsize, barrI, barrg, barcc, IdivRg, numframes, masscut);


fancytrack(basepath, FOV, featsize, maxdisp, goodenough, memory);

%%
dedrifted = dedrift(basepath, false);
%%
final = Strain_calc(dedrifted, p_size, false, 3);


%% plotting
frame_to_check = 45; % frame number
frame_to_check = frame_to_check-start_frame+1; %real frame number

%check_tracking(basepath, frame_to_check, featsize);

%truncate final array
final_truncated = {numframes};
for i = 1:numframes
    %final_truncated{i} = final(final(:,2,i)>=300 & final(:,1,i) > 100 & final(:,1,i) < 500 ,:,i);
    final_truncated{i} = final(final(:,2,i)>=300 ,:,i);
end

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'strain.gif';
for n = 1:numframes
    scatter(final_truncated{n}(:,1), -1*final_truncated{n}(:,2), 20 , final_truncated{n}(:,3),'filled');
    colorbar;
    caxis([-0.05, 0.05]);
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

%% Plot process zone for frames of interest
%final = final;  %the xx xy yy yx strain matrix
t = frame_to_check;          %the frame you want to see
strain_type = 8 ;%choose 3/4/5/6: 3 = xx; 4 = xy; 5 = yy 6 = yx
strain_cut = 0.02; %choose from below which strain value you want to exclude particles

% [circle, a, P] = calc_processZone(final,frame_to_check,strain_type,strain_cut);
% 
% figure
% plot(final_truncated(:,1,t),-1*final_truncated(:,2,t),'.');
% hold on
% plot(P(circle,1),-1*P(circle,2));
% 
% title("process Zone t="+t+", area="+a)
% xlabel("X")
% ylabel("Y")
% legend("particle", "process zone")

areas = zeros(1, numframes);
% create gif/video of all frames(if time allows)
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'processzone.gif';
v = VideoWriter('processzone', 'MPEG-4');
v.Quality = 95;
v.FrameRate = 5;
open(v);
for n = 1:numframes

    % plots the proces area at time t over all particles.
    [circle, a, P] = calc_processZone(final_truncated{n},strain_type,strain_cut);
    areas(n) = a;
    scatter(final_truncated{n}(:,1), -1*final_truncated{n}(:,2), 20 , final_truncated{n}(:,strain_type),'filled');
    colorbar;
    caxis([-0.05, 0.05]);
    hold on
    plot(P(circle,1),-1*P(circle,2));
    
    title("process Zone t="+(start_frame+n)+", area="+a)
    xlabel("X")
    ylabel("Y")
    legend("particle", "process zone")
    hold off
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 

    writeVideo(v,frame);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end
close(v);
plot(start_frame+[1:numframes], areas);