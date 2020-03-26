clear 

fname = 'data\colloids.tif';
info = imfinfo(fname);
num_images = numel(info);

images = readtiffstack(fname);
r = {num_images};

for i = 1:num_images
    disp(i)
    r{i} = feature2D(images(:, :, i),1,3,3000);
end


no = 2;

figure
imshow(images(:, :, no), [0, max(max(images(:, :, no)))]);
hold on
plot(r{2}(:,1), r{no}(:,2), '.');

remainder_x = [];
remainder_y = [];

for i = 1:num_images
    remainder_x = [remainder_x, std(mod(r{i}(:,1), 1)')];
    remainder_y = [remainder_y, std(mod(r{i}(:,2), 1)')];
end

figure
histogram(remainder_x, 10)
figure
histogram(remainder_y, 10)