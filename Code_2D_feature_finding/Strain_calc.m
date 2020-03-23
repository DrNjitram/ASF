function [final_nonzero, fig] = Strain_calc_tom(res)
% under 'Bead_tracking/res_files' find the .mat file that contains all the
% tracks. Load this.
%load('res_fov1.mat')
% rearrange the 'res' matrix to comply with Dmitry's format
MTedge = [res(:,1),res(:,2),res(:,3),res(:,7),res(:,8)];
% now go forward with his code
t1=1;
incremental=1;
%% Choose region
cutedge_p=00;
x_min=cutedge_p;x_max=max(MTedge(:,1))-cutedge_p;
y_min=00;y_max=max(MTedge(:,2));
z_min=cutedge_p;z_max=max(MTedge(:,3))-cutedge_p;
MTtrack=MTedge; %MTedge is your initial array with x,y,z coordinates, time and particle number.
MTtrack=MTtrack(MTtrack(:,1)>x_min&MTtrack(:,1)<x_max,:);
MTtrack=MTtrack(MTtrack(:,2)>y_min&MTtrack(:,2)<y_max,:);
MTtrack=MTtrack(MTtrack(:,3)>z_min&MTtrack(:,3)<z_max,:);

% Modify for 2D
MTtrack(:,3) = 1;
max_time=max(MTtrack(:,4));

%res_fig=track(MTtrack(:,1:4),6,0,max_time,3,0);
% res_fig=res_fig_c;
% %% Check track lenght distribution
% A_check=[];
% for i_n=1:max(res_fig(:,5))
%     A_check=[A_check,size(res_fig(res_fig(:,5)==i_n),1)];
% end
% figure(31);hist(A_check,[0:1:max_time])
%% Chose only the particles with full life cycle
res_t=[];
res_fig = MTtrack;
for i_n=1:max(res_fig(:,5))
    if size(res_fig(res_fig(:,5)==i_n),1)>=(max_time);%-0.25*max_time)
        res_t=[res_t;res_fig(res_fig(:,5)==i_n,:)];
    end
end


%% Choose the same particles at different times
p_size=40; %originally 40
part_n=(size(res_t,1)/max_time);
str_out=zeros(max_time,part_n,4);
D2=zeros(max_time,part_n);

for i_t=2:max_time
    t2=i_t;
    if incremental==1
        res_t1=res_t(res_t(:,4)==t2-1,:);
        res_t1(:,4:5)=[];
    else
        res_t1=res_t(res_t(:,4)==t1,:);
        res_t1(:,4:5)=[];
    end
    %     res_t1e=res_t(res_t(:,4)==t1+1,:);res_t1e(:,4:5)=[];
    %     res_t1ee=res_t(res_t(:,4)==t1+2,:);res_t1ee(:,4:5)=[];
    res_t2=res_t(res_t(:,4)==t2,:);
    res_t2(:,3:4)=[];
    %     res_t2e=res_t(res_t(:,4)==t2+1,:);res_t2e(:,4:5)=[];
    %     res_t2ee=res_t(res_t(:,4)==t2+2,:);res_t2ee(:,4:5)=[];
    %     str_aver=0;
    str_cor_t1=[res_t1(:,1:2)];
    str_cor_t2=[res_t2(:,1:2)];
    [str_out(i_t,:,:),D2(i_t,:)]=localstrain_NonAffine_mod(str_cor_t1,str_cor_t2,p_size*1.0);
    
    %average strain over neighboring particles
    
    
    
    for i_n=1:part_n
        p_d_local=sqrt((str_cor_t2(i_n,1)-str_cor_t2(:,1)).^2+(str_cor_t2(i_n,2)-str_cor_t2(:,2)).^2);%+(str_cor_t2(i_n,3)-str_cor_t2(:,3)).^2); %total distance
        c_temp=p_d_local<=p_size*2;
        str_out(i_t,i_n,:)=mean(str_out(i_t,c_temp,:),2);
        D2(i_t,i_n)=mean(D2(i_t,c_temp),2);
    end
    
   
    
end
%% TOM ADDED 12-01-17
%renumber each bead that last for the entire length of movie starting from
%1 instead of whatever the original bead# was...
% counter = 1;
% res_t(1,6) = 1;
% for i = 2:length(res_t);
%     if res_t(i,5)~=res_t(i-1,5);
%         counter = counter+1;
%         res_t(i,6) = counter;
%     else 
%         res_t(i,6) = counter;
%     end
% end

%new res_t strain component xx is column 6, xy col 7, etc...
counter = 1;
for n=1:size(str_out,2);
   for t = 1:size(str_out,1);
       if  res_t(((n-1)*size(str_out,1))+t,5) == counter;
           res_t(((n-1)*size(str_out,1))+t,6) = str_out(t,n,1); % xx component
           res_t(((n-1)*size(str_out,1))+t,7) = str_out(t,n,2); % xy component
           res_t(((n-1)*size(str_out,1))+t,8) = str_out(t,n,3); % yx component
           res_t(((n-1)*size(str_out,1))+t,9) = str_out(t,n,4); % yy component
           res_t(((n-1)*size(str_out,1))+t,10) = D2(t,n); % non-affine
           
           %display('good')
           %pause(0.1)
       else
           %display('error')
           %pause(0.1)
           counter = counter + 1; % skip bad particles
       end
   end
   counter = counter+1;
end


% rearrange to time points into a new array for xy component
final = zeros(size(str_out,2),6,size(str_out,1)); % arrange as (1) x,(2) y, (3) strain xx (4) strain xy for each time
for t = 1:size(str_out,1);
    counter = 1;
    for k=1:size(res_t,1);
        if res_t(k,4)==t;
            final(counter,1,t)=res_t(k,1); % x location
            final(counter,2,t)=res_t(k,2); % y location
            final(counter,3,t)=res_t(k,6); % xx strain
            final(counter,4,t)=res_t(k,7); % xy strain
            final(counter,5,t)=res_t(k,8); % yx strain
            final(counter,6,t)=res_t(k,9); % yy strain
            final(counter,7,t)=res_t(k,10); % non-affine
            counter = counter+1;
        end
    end
end

ypos = final(:,2,2);
ypos = abs(ypos - max(ypos));
final(:,2,2) = ypos;

% Remove all zero elements of final
[row_xy_strain,~]=find(final(:,3,2));
%
final_nonzero(:,1,2) = final(row_xy_strain,1,2);
final_nonzero(:,2,2) = final(row_xy_strain,2,2);
final_nonzero(:,3,2) = final(row_xy_strain,3,2);
final_nonzero(:,4,2) = final(row_xy_strain,4,2);
final_nonzero(:,5,2) = final(row_xy_strain,5,2);
final_nonzero(:,6,2) = final(row_xy_strain,6,2);
final_nonzero(:,7,2) = final(row_xy_strain,7,2);
%

