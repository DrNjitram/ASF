function [circle, a, P] = calc_processZone(final,strain_type,strain_cut)
%input is the final matrix from the strain_calc.m function
%output is 'beadID'x2 matrix with the positions of beads over a certain 
    if strain_cut == -1
       strain_cut = min(final(:,strain_type)) + 0.8 * (max(final(:,strain_type)) - min(final(:,strain_type)));
    end
    P = final(final(:,strain_type)>=strain_cut,1:2); % ID = :, strain_type = 3-6, t =
    if length(P) > 2
        [circle, a] = convhull(P); % a = area of convex hull
    else
        circle = [];
        a = 0;
    end
end