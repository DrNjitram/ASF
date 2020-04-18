function a = calc_processZone(final,t,strain_type,strain_cut)
%input is the final matrix from the strain_calc.m function
%output is 'beadID'x2 matrix with the positions of beads over a certain 
P = final(final(:,strain_type,t)>=strain_cut,1:2,1); % ID = :, strain_type = 3-6, t =
[circle, a] = convhull(P); % a = area of convex hull

% plots the proces area at time t over all particles.
figure
plot(final(:,1,t),final(:,2,t),'.');
hold on
plot(P(circle,1),P(circle,2));
title("process Zone t="+t+", area="+a)
xlabel("X")
ylabel("Y")
legend("particle", "process zone")