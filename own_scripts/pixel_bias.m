function [] = pixel_bias(res)
    %% Calculate standard deviation of remainders mod 1

    remainder_x = mod(res(:,1), 1)';
    remainder_y = mod(res(:,2), 1)';

    %% Imaging
    % Display histograms of remainders
    figure
    histogram(remainder_x, 10); hold on;
    histogram(remainder_y, 10);


end