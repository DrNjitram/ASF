% convhull(P);
% input is een xy matrix.
% De waarden in de xy matrix moeten voldoen aan strain>iets
% ergo
% we moeten waarden uit 'final' pakken die voldoen aan final(ID,kolom 3-6,t)>iets:

%% just a setup for now
% for loop to calculate process zone for every time
for t = 1:max(final(:,:,3)) %misschien werkt max(final,3) wel.
    % selects correct x and y values with a corresponding strain higher
    % than something
    P{t} = final(final(ID,strain_type,t)>=iets,:); % ID = :, strain_type = 3-6, t = loop 
end
% calculates the convex hul and its surface area for time t
[circle, a] = convhull(P(t)); % a = area of convex hull

% plots the proces area at time t over all particles.
plot(final(:,1,t),final(:,2,t),'*');
hold on
plot(final(P,1),final(P,2));
