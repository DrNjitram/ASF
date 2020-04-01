%%create 2 datasets for different frames
M1 = randn(100, 8);
M2 = [M1(:,1), M1(:,2)];
M0 = [ones(100,1),ones(100,1)*2];
M3 = M2+M0;

MT1 = [M2(:,:),ones(100,1)]; %x and y, frame 1
MT2 = [M3(:,:),ones(100,1)*2]; %x+1 and y+1, frame 2


%%Should be done for every frame
F1 = [MT1(:,1),MT1(:,2)];
F2 = [MT2(:,1),MT2(:,2)];
delta = F2-F1;
Mx = [median(delta(:,1))];
My = [median(delta(:,2))];
MedianVector = [Mx*ones(100,1) ,My*ones(100,1), zeros(100,1)];
MT2 = [MT2(:,1)-MedianVector(:,1),MT2(:,2)-MedianVector(:,2),MT2(:,3)];