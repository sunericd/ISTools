function [centers] = CellCount(fname, tp, N)
    % fname = filename of image
    % tp = timepoint; corresponds to layer of image.
    % N = sensitivity
    
    % Finds and reads the image
    img = imread(fname, tp);

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
    centers = regionprops(erod, 'Centroid');
end