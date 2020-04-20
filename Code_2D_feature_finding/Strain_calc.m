function [final] = Strain_calc(res, p_size, fig, type)
% Under 'Bead_tracking/res_files' find the .mat file that contains all the
% tracks. Load this data first.
% This code was written for 3D data and has been adapted to your 2D data
% set
% Most of this code is about rearranging data to consider only the
% neightbors of the target particle from which to calculate the local
% strain
fprintf('Calculating Strain:\n');
% rearrange the 'res' matrix to comply with strain calculation's format
MTedge = [res(:,1),res(:,2),res(:,3),res(:,6),res(:,8)];
% now go forward with this code
t1=1;
frame_time = 4.2; %change this as necessary
incremental=1;

%% Choose region of interest
cutedge_p=0; %choose a number of pixels to ignore from the edges of the image; now set to include the entire image
x_min=cutedge_p;
x_max=max(MTedge(:,1))-cutedge_p;

y_min=0;
y_max=max(MTedge(:,2));

z_min=cutedge_p; %There is no z value in 2D data, this will be fixed a few lines later in the code.
%z_max=max(MTedge(:,3))-cutedge_p;
z_max=max(MTedge(:,1))-cutedge_p; %set max to largest dimension in y, again this does not matter since we have no z-data

MTtrack=MTedge; %MTedge is your initial array with x,y,z coordinates, time and particle number.
MTtrack=MTtrack(MTtrack(:,1)>x_min&MTtrack(:,1)<x_max,:);
MTtrack=MTtrack(MTtrack(:,2)>y_min&MTtrack(:,2)<y_max,:);

% Modify for 2D
MTtrack(:,3) = 1;
%MTtrack(:,4) = res(:, 6); %change the real time in second into frame number as this makes the next 50 lines easier to code in whole integers; use 'ceil' to round to next highest integer
max_time=max(MTtrack(:,4)); %Note: MTtrack(:,4) == res(:,7)

%% Chose only the particles with full life cycle in the chosen region
% If you set 'goodenough' equal to the full number of frames, this code is
% not important, but run it still as the variable names change. 
res_t=[];
res_fig = MTtrack;
fprintf('Selecting Features ');
linelength = 0;
targ = max(res_fig(:,5)); 
for i_n=1:targ %Note: res_fig(:,5) are BeadID, so i_n goes from 1 to max(BeadIDs)
    fprintf(repmat('\b',1,linelength));
    linelength = fprintf('%d of %d', i_n, targ );
    if size(res_fig(res_fig(:,5)==i_n),1)>=(max_time) %use if to select only particles with frames = max_frame
        res_t=[res_t;res_fig(res_fig(:,5)==i_n,:)]; %make that selection, ignore all other BeadID; store as new variable "res_t". the first part of this line that say [res(t)...] causes it to grow each loop
    end
end
fprintf('\n');

%% Choose the same particles at different times and calculat the strain
% this takes some time as there are many frames and many particles
%p_size=10; %average distance between neighboring particles in pixels; choose a value much greater than diameter; this is an input to the strain code

part_n=(size(res_t,1)/max_time); %the size of res_t / max_time is the number of total particles
str_out=zeros(max_time,part_n,4); %create an empty array for strain (i.e. str_out); the 4 here refers to a 4-component strain tensor [xx, xy, yy, yx]

D2=zeros(max_time,part_n); %create an empty array for non-affinitiey (i.e. D2)

fprintf('Processing frame ')
linelength = 0;
for i_t=2:max_time %start at frame = 2
    fprintf(repmat('\b',1,linelength));
    linelength = fprintf('%d of %d', i_t, max_time);
    
    t2=i_t; %define the 'second' frame; this changes during the for-loop
    
    res_t1=res_t(res_t(:,4)==t2-1,:); %choose the res_t values for the frame "t2-1"
    res_t1(:,4:5)=[]; %get rid of not needed variables
    
    res_t2=res_t(res_t(:,4)==t2,:); %choose the res_t values for the frame "t2"
    res_t2(:,3:4)=[]; %get rid of not needed variables

    
    str_cor_t1=[res_t1(:,1:2)]; %basically just rename the x,y,z position of the particles to this variable name for the first frame
    str_cor_t2=[res_t2(:,1:2)]; %basically just rename the x,y,z position of the particles to this variable name for the second frame
    [str_out(i_t,:,:),D2(i_t,:)]=localstrain_2D(str_cor_t1,str_cor_t2,p_size); %This code calculate each particle's strain for the two frames
   
    
    
    %average strain over neighboring particles
       
    for i_n=1:part_n
        p_d_local=sqrt((str_cor_t2(i_n,1)-str_cor_t2(:,1)).^2+(str_cor_t2(i_n,2)-str_cor_t2(:,2)).^2);%total distance
        c_temp=p_d_local<=p_size*2; %determines neighbor particles within p_size
        str_out(i_t,i_n,:)=mean(str_out(i_t,c_temp,:),2); %takes the mean values of strain for those neighbor particles
        D2(i_t,i_n)=mean(D2(i_t,c_temp),2); %takes the mean values of non-affinity (D2) for those neighbor particles
    end
    
   
    
end
fprintf('\n\n');
%% TOM ADDED 31-03-2020
% Now add the strain components to the x, y, positions found in variable
% 'res_t'...

% Strain in 2D is a 4 component vector
%new res_t strain component xx is column 6, xy col 7, etc...
counter = 1;
for n=1:size(str_out,2) %loop over number of particles
   for t = 1:size(str_out,1) %loop over frame number 
       if  res_t(((n-1)*size(str_out,1))+t,5) == counter
           res_t(((n-1)*size(str_out,1))+t,6) = str_out(t,n,1); % xx component
           res_t(((n-1)*size(str_out,1))+t,7) = str_out(t,n,2); % xy component
           res_t(((n-1)*size(str_out,1))+t,8) = str_out(t,n,3); % yx component
           res_t(((n-1)*size(str_out,1))+t,9) = str_out(t,n,4); % yy component
           res_t(((n-1)*size(str_out,1))+t,10) = D2(t,n); % non-affine
           

       else
           counter = counter + 1; % skip particle IDs that are missing from res_t(:,5)
       end
   end
   counter = counter+1;
end

%% Reduce data complexity
% rearrange to time points into a new array of xy component, 4 strains,
% and D2 called 'final'
final = zeros(size(str_out,2),6,size(str_out,1)); % arrange as (1) x,(2) y, (3) strain xx (4) strain xy for each time
for t = 1:size(str_out,1)
    counter = 1;
    for k=1:size(res_t,1)
        if res_t(k,4)==t
            final(counter,1,t)=res_t(k,1); % x location
            final(counter,2,t)=res_t(k,2); % y location
            final(counter,3,t)=res_t(k,6); % xx strain
            final(counter,4,t)=res_t(k,7); % xy strain
            final(counter,5,t)=res_t(k,8); % yx strain
            final(counter,6,t)=res_t(k,9); % yy strain
            final(counter,7,t)=res_t(k,10); % non-affine
            final(counter,8,t)=res_t(k,6)+res_t(k,7)+res_t(k,8)+res_t(k,9); % all combined
            final(counter,9,t)=res_t(k,7)+res_t(k,8); %diagonals
            counter = counter+1;
        end
    end
end


%% Time to make a plots...
%frame = 78;
%scatter(res(res(:,7)==frame, 1), res(res(:,7)==frame, 2), 20, zeros(length(res(res(:,7)==frame, 1)), 1));
%scatter(final(:,1,frame), final(:,2,frame), 20 , final(:,4,frame),'filled'); hold on; %this will make a scatter plot of x-, y- positions, with each point having a size = 20, and color of th xy strain component for only one frame
if fig
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'strain.gif';
    for n = 1:max_time


        scatter(final(:,1,n), -1*final(:,2,n), 20 , final(:,type,n),'filled');
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
end
