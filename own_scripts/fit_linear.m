function [remainder, values] = fit_linear(dataset, x, limit)
    if nargin < 4; limit = 3; end
    
    guess_ax = (max(dataset)-min(dataset))/length(x);
    guess_bx = dataset(fix(length(x)/2))-fix(length(x)/2)*guess_ax;

    drifta = @(a, x) a*x + guess_bx;
    [~, a] = find_min(drifta, x, guess_ax, [0, 10], dataset);

    driftb = @(b, x) a*x + b;
    [remainder, b] = find_min(driftb, x, guess_bx, [0, 10], dataset);
    drifta = @(a, x) a*x + b;
    
    for i = 0:limit
        [~, a] = find_min(drifta, x, a, [0, 10], dataset);
        [remainder, b] = find_min(driftb, x, b, [0, 10], dataset);
    end
    values = [a, b];
end