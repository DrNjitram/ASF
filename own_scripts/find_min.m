function [remainder, value] = find_min(fit_func, x, guess, range, target,  steps, recursive_depth_max, l_remainder, fig, n)
    if nargin < 6; steps = 1000; end    
    if nargin < 7; recursive_depth_max = 3; end
    if nargin < 8; l_remainder = inf; end
    if nargin < 9; fig = false; end
    if nargin < 10; n = 1; end
    
    rs = zeros(steps, 1);
    vs = zeros(steps, 1);
    
    if fig
        figure;
        axis([range(1)*guess range(2)*guess 0 1]);
    end
    
    step = (range(2) * guess - range(1) * guess)/steps;
    
    for iter = 0:steps
        a = guess * range(1) + step * iter;
        S = fit_func(a, x);
        remainder = sum(sqrt(mean((target - S).^2)));

        rs(iter+1) = remainder;
        vs(iter+1) = a;
        
        if fig
            plot(a, remainder, 'b.'); hold on;
        end
    end
    
    [remainder, index] = min(rs);
    value = vs(index);
    
    if remainder ~= l_remainder && recursive_depth_max > 0 
         [remainder, value] = find_min(fit_func, x, value, [1-10^-n, 1+10^-n], target,  steps, recursive_depth_max - 1, remainder, fig, n+1);
    end
end