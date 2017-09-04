function [cell_score] = CellLite(dir_path1, dir_path2, extension1, extension2, N)

% dirpath_1 and dirpath_2 are directory paths ending with '/' and have
% images ordered by comparison (have same named components)

intensity_score = [];
cell_ratio = [];
cell_score = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 = dapi, 2 = tunel --> does not have to be a dapi or tunel image
dapi_imagenames = dir(fullfile(dir_path1, extension1));
tunel_imagenames = dir(fullfile(dir_path2, extension2));

length(dapi_imagenames)


for i = 1:length(dapi_imagenames)
    
    %experiments(i) = dapi_imagenames(i).name
    
    % INITIAL : PScG_100_2_1_DAPI/FITC

    image_dapi = imread(fullfile(dir_path1, dapi_imagenames(i).name));
    image_tunel = imread(fullfile(dir_path2, tunel_imagenames(i).name));

    bw_image_dapi = image_dapi(:,:,3);
    bw_image_tunel = rgb2gray(image_tunel);
    
    % Thresholding
    move = 3;
    up_d = imtranslate(bw_image_dapi, [0, -move]);
    left_d = imtranslate(bw_image_dapi, [-move, 0]);
    up_t = imtranslate(bw_image_tunel, [0, -move]);
    left_t = imtranslate(bw_image_tunel, [-move, 0]);

    up_subtract_d = im2double(bw_image_dapi) - im2double(up_d);
    left_subtract_d = im2double(bw_image_dapi) - im2double(left_d);
    up_subtract_t = im2double(bw_image_tunel) - im2double(up_t);
    left_subtract_t = im2double(bw_image_tunel) - im2double(left_t);

    roughness_d = up_subtract_d.^2 + left_subtract_d.^2;
    roughness_d = roughness_d/max(max(roughness_d));
    roughness_t = up_subtract_t.^2 + left_subtract_t.^2;
    roughness_t = roughness_t/max(max(roughness_t));

    thresh_image_dapi = im2bw(sqrt(roughness_d(1:(size(roughness_d, 1)-move), (move+1):(size(roughness_d, 2)))), 0.01);
    thresh_image_dapi = im2bw(roughness_d, 0.008);
    thresh_image_tunel = im2bw(sqrt(roughness_t(1:(size(roughness_t, 1)-move), (move+1):(size(roughness_t, 2)))), 0.01);
    thresh_image_tunel = im2bw(roughness_t, 0.008);
    
    % Cleaning (Fill and Erode)
    se22 = strel('disk', round(3*N));
    se23 = strel('disk', round(8*N));

    thresh_image_dapi = imdilate(thresh_image_dapi, se22);
    thresh_image_dapi = imfill(thresh_image_dapi, 'holes');
    thresh_image_dapi = imerode(thresh_image_dapi, se23);
    
    thresh_image_tunel = imdilate(thresh_image_tunel, se22);
    thresh_image_tunel = imfill(thresh_image_tunel, 'holes');
    thresh_image_tunel = imerode(thresh_image_tunel, se23);

    %figure; imshow(thresh_image_dapi)
    %figure; imshow(thresh_image_tunel)
    
    % Labeling
    label_image_dapi = bwlabel(thresh_image_dapi);
    label_image_tunel = bwlabel(thresh_image_tunel);

    blob_num_dapi = max(max(label_image_dapi)); % gives you the number of blobs
    blob_num_tunel = max(max(label_image_tunel)); % gives you the number of blobs


    blob_sizes_dapi = zeros(blob_num_dapi, 1); % list of blob sizes
    blob_sizes_tunel = zeros(blob_num_tunel, 1); % list of blob sizes


    for j = 1:blob_num_dapi; % populates blob_sizes
        frq = length(find(label_image_dapi==j));
        blob_sizes_dapi(j) = frq;
    end

    for k = 1:blob_num_tunel; % populates blob_sizes
        frq = length(find(label_image_tunel==k));
        blob_sizes_tunel(k) = frq;
    end

    max_thresh_dapi = 0.01*max(blob_sizes_dapi); % threshold for what is big enough for a cell (adjusted for DAPI)
    max_thresh_tunel = 0.04*max(blob_sizes_tunel); % threshold for what is big enough for a cell (adjusted for Tunel)


    cell_sizes_dapi = []; 
    cell_sizes_tunel = [];

    % DAPI
    cell_centers_y_dapi = [];
    cell_centers_x_dapi = [];

    for l = 1:length(blob_sizes_dapi);
        if blob_sizes_dapi(l) > max_thresh_dapi;
            [x, y] = find(label_image_dapi==l);
            cell_centers_x_dapi(length(cell_centers_x_dapi)+1) = round(mean(x));
            cell_centers_y_dapi(length(cell_centers_y_dapi)+1) = round(mean(y));
            cell_sizes_dapi(length(cell_sizes_dapi)+1) = blob_sizes_dapi(l);
        end
    end

    % Tunel

    cell_centers_y_tunel = [];
    cell_centers_x_tunel = [];

    for m = 1:length(blob_sizes_tunel);
        if blob_sizes_tunel(m) > max_thresh_tunel;
            [x, y] = find(label_image_tunel==m);
            cell_centers_x_tunel(length(cell_centers_x_tunel)+1) = round(mean(x));
            cell_centers_y_tunel(length(cell_centers_y_tunel)+1) = round(mean(y));
            cell_sizes_tunel(length(cell_sizes_tunel)+1) = blob_sizes_tunel(m);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%

    % NUMBER CELLS
    cell_num_dapi = length(cell_sizes_dapi);
    cell_num_tunel = length(cell_sizes_tunel);

    cell_diff = cell_num_dapi - cell_num_tunel;
    cell_ratio(i) = cell_num_tunel/cell_num_dapi;

    % INTENSITY (proportional to sizes)
    intensity_dapi = sum(cell_sizes_dapi);
    intensity_tunel = sum(cell_sizes_tunel);

    intensity_score(i) = intensity_tunel/intensity_dapi;
    
    % CELL SCORE (Intensity/Num)
    cell_score(i) = (intensity_tunel/intensity_dapi)/(cell_num_tunel/cell_num_dapi);

end
end
