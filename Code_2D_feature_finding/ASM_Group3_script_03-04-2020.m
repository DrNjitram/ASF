
load('locating_variables.mat')
featsize = 3;
barrI = 1000;
Imin = 100;
mpretrack(basepath,fovn,featsize,barrI,barrRg,barcc,IdivRg,numframes,masscut,Imin,2)


load(['.\Feature_finding\MT_1_Feat_Size_',num2str(floor(featsize)),'.mat']) %concatenate strings using featsize variable
frames = find(MT(:,6)==2); %find where (which element) where the second frame starts

%%
%check visually if particle locations match the first frame
first_image = imread('.\fov1\fov1_0001.tif');
imagesc(first_image); hold on;
scatter(MT(1:(frames(1)-1),1),MT(1:(frames(1)-1),2),'*k'); 




%%
%check pixel biasing
xbias = mod(MT(:,1),1); %find using modulus, find the remainder
ybias = mod(MT(:,2),1); 

histogram(xbias,10,'FaceColor','r'); hold on;
histogram(ybias,10,'FaceColor','m');
%%
% tracking 
fancytrack( basepath, fovn, featsize, maxdisp, goodenough, memory ) %all inputs required


%%
% Plot 50 of those particles
for m=1:50;
    plot(res((1+(79*(m-1))):(79*m),1),res((1+(79*(m-1))):(79*m),2)); hold on;
end
%%
% dedrifting code
% explained poorly in description document; take the average position of
% all particles at t=2, and subtract the previous average position at t=1;
% this will give the displacement due to drift. Subtract this displacement
% from each particle position at t=2. 

%%
% calculate local strain for each particle