function [num_cells] = CellCountInitial(directory_path, extension, N)
    % directory_path = string path to folder (ends with \ or /)
    % extension = file format extension (e.g. *.tif, *.jpg)
    % N = sensitivity

    % Finds the files in the current directory that share the
    % directory and extension.
    imagenames = dir(fullfile(directory_path, extension));

    
    % Finds and reads the image
    filename = imagenames(1).name;
    img = imread(fullfile(directory_path, filename));

    % Turns RGB to BW
    bw_image = rgb2gray(img);

    % CELL SHAPE
    % Thresholding
    wb_image = imcomplement(bw_image); % reverses B and W

    % Directional Derivative Thresholding
    move = 3;
    up = imtranslate(wb_image, [0, -move]);
    left = imtranslate(wb_image, [-move, 0]);

    up_subtract = im2double(wb_image) - im2double(up);
    left_subtract = im2double(wb_image) - im2double(left);

    roughness = up_subtract.^2 + left_subtract.^2;
    roughness = roughness/max(max(roughness));

    derivate = im2bw(sqrt(roughness(1:(size(roughness, 1)-move), (move+1):(size(roughness, 2)))), 0.01);
    derivate = im2bw(roughness, 0.008);

    % Dilating and Filling Eroding

    se22 = strel('disk', round(3*N));
    se23 = strel('disk', round(8*N));

    dilate = imdilate(derivate, se22);
    dilate = imfill(dilate, 'holes');
    erod = imerode(dilate, se23);

    % Finding the centers
    centers = regionprops(erod, 'Centroid')
    num_cells = length(centers);
    
    % Plotting centeres on initial image...
    figure; imshow(img,[]); colormap(jet); hold on;
    
    for i = 1:length(centers)
        plot(centers(i).Centroid(1), centers(i).Centroid(2),'ro');
    end
        
   
end

    
