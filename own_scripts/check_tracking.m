function [] = check_tracking(basepath, frame_to_check, featsize)
    load([basepath '\Feature_finding\MT_1_Feat_Size_',num2str(floor(featsize)),'.mat']) %concatenate strings using featsize variable
    load( [basepath 'Bead_tracking/res_files/res_fov1.mat'] ); %load relevant data file

    frames = find(MT(:,6)==frame_to_check); %find where (which element) where the second frame starts
    first_image = imread([basepath '\fov1\fov1_' num2str(frame_to_check,'%04i') '.tif']);
    imagesc(first_image); hold on;
    scatter(MT((frames),1), MT((frames),2),'*k'); 

    %% Calculate standard deviation of remainders mod 1

    remainder_x = mod(res(:,1), 1)';
    remainder_y = mod(res(:,2), 1)';

    %% Imaging
    % Display histograms of remainders
    figure
    histogram(remainder_x, 10); hold on;
    histogram(remainder_y, 10);


end