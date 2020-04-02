function subimages = readtiffstack(filename)
    num_images = numel(imfinfo(filename));
    
    subimages = {num_images};
    
    for i = 1:num_images
        subimages{i} = imread(filename, i);
    end
end