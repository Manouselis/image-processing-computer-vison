function [disparityMap, reliability_map] = calculateDisparityMap(img_pair,...
    disparityRange, reference)
% Convert images to greyscale
img_1_grey = rgb2gray(img_pair{1});
img_2_grey = rgb2gray(img_pair{2});

% If reference = 2 (disparity map to be calculate wrt 2nd image), then:
if reference == 2
    % disparity range has to be negative
    disparityRange = -disparityRange;
    disparityRange([1 2]) = disparityRange([2 1]);

    % the order of images has to be swaped
    img_1_grey = rgb2gray(img_pair{2});
    img_2_grey = rgb2gray(img_pair{1});
end

% Calculate disparity map
disparityMap = disparitySGM(img_1_grey, img_2_grey,...
    'DisparityRange',disparityRange,'UniquenessThreshold',0);

% If disparity map was calcualted wrt 2nd image, then disparities are
% negative. They are to be multiplied by -1 to convert them to positive:
if reference == 2
    disparityMap = -disparityMap;
end
reliability_map = ~isnan(disparityMap);
reliability_map(disparityMap<0) = 0;
end