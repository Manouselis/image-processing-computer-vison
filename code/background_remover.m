function face_mask = background_remover(im)
%% Step 1 - Saturation thresholding for background detection
% The subjects have been photographed in a room where the color of the walls
% are close to the color of their skin
% For better color segmentation firstly we convert the RGB image to an HSV
% image and threshold the saturation level
% Convert from RGB to HSV color space
im_hsv = rgb2hsv(im);
% Split the channels
s = im_hsv(:, :, 2);
% Get number of bins and values for them
[N, edges] = histcounts(s);
idx = find(N == max(N));

% figure; histogram(s);

% Get the first local minimum and it's index
lmin = islocalmin(N(:, idx:size(N, 2)));
lmin_idx = find(lmin == 1);
% Don't get the first value, have some offset
lmin_idx = lmin_idx(5);

% Get the threshold value
threshold_value = edges(:, idx + lmin_idx);

% Create mask based on chosen histogram thresholds
s_mask = (im_hsv(:,:,1) >= 0 ) & (im_hsv(:,:,1) <= 1) & ...
    (im_hsv(:,:,2) >= threshold_value ) & (im_hsv(:,:,2) <= 1) & ...
    (im_hsv(:,:,3) >= 0 ) & (im_hsv(:,:,3) <= 1);

% Initialize output masked image based on input image
im_masked = im;

% Set background pixels where BW is false to zero
im_masked(repmat(~s_mask,[1 1 3])) = 0;
% figure; imshow(im_masked);

%% Step 2 - Color based segmentation using K-means to find the person
% The first color based segmentation is used to find the person, meaning
% face + outfit
% Convert the input image into L*a*b* space color
im_lab = rgb2lab(im);

% Extract a*b*
ab = im_lab(:, :, 2:3);
ab = im2single(ab);

% Get the pixel labels for 3 colors (background, outfit and face)
pixel_labels = imsegkmeans(ab, 3);

% The pixel label with the second highest ocurances should be the face, the
% third one is the outfit
counted_pixel_labels = [sum((pixel_labels(:) == 1));
    sum((pixel_labels(:) == 2));
    sum((pixel_labels(:) == 3))];

sorted_pixel_labels = sort(counted_pixel_labels,'descend');

% Extract the face and outfit pixels
face_pixels = find(counted_pixel_labels == sorted_pixel_labels(2));
outfit_pixels = find(counted_pixel_labels == sorted_pixel_labels(3));

face_mask = (pixel_labels == face_pixels);
outfit_mask = (pixel_labels == outfit_pixels);

% Create a binary mask to for the person
face_mask = face_mask | outfit_mask;

% figure; imshow(labeloverlay(im,pixel_labels));

%% Step 3 - Perform morphological operations to get the mask of the person
% Extract the greatest area from the opened face mask which is the face
face_mask = bwareafilt(face_mask, 1);
% figure; imshow(face_mask);

% Fill in the holes and gaps in the mask
face_mask = imfill(face_mask,'holes');
% figure; imshow(face_mask);

%
face = im;
face(repmat(face_mask,[1,1,3])==0) = 0;
% result = labeloverlay(im,face_mask);
% figure; imshow(result);
im(repmat(face_mask,[1,1,3])==0)=0;
% figure; imshow(im);

%% Step 4 - Color based segmentation using K-means for face extraction
% Now that the blackground is completly black, the face and outifit of the
% subjects can be separted more easily
% Convert the input image into L*a*b* space color
im_lab = rgb2lab(im);

% Extract a*b*
ab = im_lab(:, :, 2:3);
ab = im2single(ab);

% Get the pixel labels for 3 colors (background, outfit and face)
pixel_labels = imsegkmeans(ab, 3);

% The pixel label with the second highest ocurances should be the face
counted_pixel_labels = [sum((pixel_labels(:) == 1));
    sum((pixel_labels(:) == 2));
    sum((pixel_labels(:) == 3))];

sorted_pixel_labels = sort(counted_pixel_labels,'descend');

% Extract the face pixels
face_pixels = find(counted_pixel_labels == sorted_pixel_labels(2));
face_mask = (pixel_labels == face_pixels);

% figure; imshow(labeloverlay(im,pixel_labels));

%% Step 5 - Perform morphological operations to get the mask of subject's face
% Extract the greatest area from the opened face mask which is the face
face_mask = bwareafilt(face_mask, 1);
% figure; imshow(face_mask);

se = strel('disk', 5);
face_mask = imerode(face_mask, se);
% figure; imshow(face_mask);

% Fill in the holes and gaps in the mask
face_mask = imfill(face_mask,'holes');
% figure; imshow(face_mask);

% Dilate the image with horizontal and vertical lines to ensure
% connectivity
se = strel('line', 10, 90);
se_horizontal = strel('line', 10, 0);
face_mask = imdilate(face_mask, [se se_horizontal]);
% figure; imshow(face_mask);

% Just to be sure
face_mask = imfill(face_mask,'holes');

% imwrite(face_mask, 'face_mask.jpg');

% Get the convexhull of the of the dilated mask
% face_mask = bwconvhull(face_mask, 'objects', 8);
% figure; imshow(face_mask);
end