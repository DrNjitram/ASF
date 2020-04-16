function res_dedrifted = dedrift(fig)
    %% Driftcorrection

    % INPUT
    % 
    % matrix called res, originating from 'fancytrack.m' ADD MORE DETAILS
    % CONCERNING MATRIX FORMAT

    % OUTPUT
    % 
    % 

    load('res_fov1.mat') %load relevant data file
    res_yes = [res(:,1),res(:,2),res(:,6),res(:,8)]; %remove irrelevant values

    %% drift correction based on average positions
    % INPUT: [res_yes]
    % OUTPUT: []

    sortFOV = sortrows(res_yes, 3); %sorts the matrix to allow easy selection of frame i
    corMatrix = zeros(length(res_yes), 2); %create array to hold subtracted averages

    frame_number = max(res_yes(:,3)); %gives final frame for defining for loop below

    meanX = zeros(1, frame_number);
    meanY = zeros(1, frame_number);
    corMeanX = zeros(1, frame_number);
    corMeanY = zeros(1, frame_number);

    for i = 1:frame_number
        FOV = find(sortFOV(:,3) == i); %gives vector containing every row belonging to frame i
        maxRow = max(FOV);
        minRow = min(FOV);
        X = sortFOV(minRow:maxRow,1); %extracts x-values belonging to frame i
        Y = sortFOV(minRow:maxRow,2); %same as above for y-values

        meanX(i) = median(X); %calculates average x-position
        meanY(i) = median(Y); %calculates average Y-position

        if i == 1
            corMeanX(1) = 0;
            corMeanY(1) = 0;
        else
            corMeanX(i) = meanX(i)-meanX(1); %subtracts average position of the frame before from the average position of a frame.
            corMeanY(i) = meanY(i)-meanY(1); %same as above but for Y
        end

        corMatrix(minRow:maxRow,1) = corMeanX(i); %puts the x-values we use for correction at the correct place
        corMatrix(minRow:maxRow,2) = corMeanY(i); %same as above for y-values
    end

    corSortFOV = [sortFOV(:,1:2)- corMatrix(:,1:2),sortFOV(:,3:4)]; %subtracts the correction value from the original positions to give the drift corrected positions
    corRes_yes = sortrows(corSortFOV,4); % sorted by beadID again

    %% Fit linear equations

    %form drift = a*frame_no + b
    frames = 1:frame_number;

    [~, values_a] = fit_linear(corMeanX, frames);
    x_fit = values_a(1) * frames + values_a(2);
    x_fit_func = @(x) values_a(1) * x + values_a(2);

    [~, values_b] = fit_linear(corMeanY, frames);
    y_fit = values_b(1) * frames + values_b(2);
    y_fit_func = @(x) values_b(1) * x + values_b(2);

    %% Apply linear equations
    res_correction = [x_fit_func(res(:, 6)), y_fit_func(res(:, 6))];
    res_dedrifted = res;
    res_dedrifted(:, 1:2) = res_dedrifted(:, 1:2) - res_correction;
    
    
    %% Next step: Show the movement plots of particles (Tom's code)
    % Plots 50 particles as a function of frame for both pre-driftcorrection
    % (figure 1) and post-driftcorrection (figure 2)
    if fig
        figure
        subplot(1, 2, 1);
        plot(frames, corMeanX); hold on;
        plot(frames, corMeanY); hold on;
        plot(frames, x_fit); hold on;
        plot(frames, y_fit); 
        xlabel("Frame Number");
        ylabel("Drift (px)");
        legend("X drift", "Y drift", "X fit", "Y fit");
        subplot(1, 2, 2);
        plot(meanX, meanY); hold on;
        plot(meanX(1) + x_fit,meanY(1) + y_fit);
        xlabel("X Drift (px)");
        ylabel("Y Drift (px)");
        legend("XY drift", "fit");

        figure 
        %pre-correction
        subplot(1, 2, 1);
        ps = 1;% particles to plot
        pe = 50;
        for m=ps:pe
            FOV = find(res(:,8) == m);
            plot(res(min(FOV):max(FOV), 1), res(min(FOV):max(FOV), 2)); hold on;
        end
        xlabel("X");
        ylabel("Y");
        title("Raw Position Data");

        %post-correction
        subplot(1, 2, 2);
        for m=ps:pe
            FOV = find(res(:,8) == m);
            plot(res_dedrifted(min(FOV):max(FOV), 1), res_dedrifted(min(FOV):max(FOV), 2)); hold on;
        end

        xlabel("X");
        ylabel("Y");
        title("Dedrifted Position Data");
    end