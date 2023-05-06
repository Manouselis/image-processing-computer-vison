function [D_map_out, reliability_map] = enhanceDmap(D_map, reliability_map,...
    img_pair, disparityRange, mask,reference)
D_map_out = D_map;
%% First enhancement: Calculate second disparity map
% To filter out points that have different disparities in both maps

if reference == 1
    % If the provided D_map is wrt first image, then calculate the
    % D_map wrt the second image
    D_map2 = calculateDisparityMap(img_pair, disparityRange,2);
    % For a point in disparity map1, calculate its disparity in disparity map2
    corresponding_disparities = D_map;
    for y=1:size(D_map,1)
        for x=1:size(D_map,2)
            if floor(x-D_map(y,x)) <= 0 || isnan(D_map(y,x))
                corresponding_disparities(y,x) = -realmax;
                continue
            end
            % notice the -ve sign in x-D_map(y,x)
            corresponding_disparities(y,x) = D_map2(y,floor(x-D_map(y,x)));
        end
    end

else % reference == 2
    % If the provided D_map is wrt second image, then calculate the
    % D_map wrt the first image
    D_map2 = calculateDisparityMap(img_pair, disparityRange,1);
    % For a point in disparity map1, calculate its disparity in disparity map2
    corresponding_disparities = D_map;
    for y=1:size(D_map,1)
        for x=1:size(D_map,2)
            if floor(x-D_map(y,x)) <= 0 || isnan(D_map(y,x))
                corresponding_disparities(y,x) = -realmax;
                continue
            end
            % notice the +ve sign in x+D_map(y,x)
            corresponding_disparities(y,x) = D_map2(y,floor(x+D_map(y,x)));
        end
    end
end
% Mark points of same disparities in both maps as reliable
reliability_map(abs(D_map - corresponding_disparities)<1) = 1;

%% Second enhancement: Any point that lie on the background is unreliable
reliability_map = reliability_map.*mask{reference};

%% Third enhancement: Remove holes in the face
reliability_map = imfill(reliability_map,'holes');
D_map_out(~reliability_map) = -realmax;

%% Fourth enhancement: Using Median Filtering
D_map_out = medfilt2(D_map_out, [40 40],'symmetric');
% Update reliability_map:
reliability_map(D_map_out==-realmax) = 0;
reliability_map(D_map_out==-Inf) = 0;
reliability_map(D_map_out==Inf) = 0;
reliability_map(isnan(D_map_out)) = 0;
D_map_out(~reliability_map) = NaN;
end