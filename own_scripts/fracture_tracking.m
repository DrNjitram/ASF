[image,map] = imread('data\colloids.tif', 1);
[r] = feature2D(image,1,3,3000);

figure
imshow(image, [0, max(max(image))]);
hold on
plot(r(:,1), r(:,2), '.');

figure
histogram(mod(r(:,1), 1), 10)
figure
histogram(mod(r(:,2), 1), 10)