function [out] = calc_processZone(t,strain_type,strain_value)
%input is the final matrix from the strain_calc.m function
%output is 'beadID'x2 matrix with the positions of beads over a certain 
P = final(final(:,strain_type,t)>=strain_value,1:2,1); % ID = :, strain_type = 3-6, t =
[circle, a] = convhull(P); % a = area of convex hull

% plots the proces area at time t over all particles.
plot(final(:,1,t),final(:,2,t),'*');
hold on
plot(P(circle,1),P(circle,2));
out = P;

%post pictures