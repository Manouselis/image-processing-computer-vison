clear variables
close all
clc
%% Preparing the environment
% Get the images of the subject.
subject = 4; % Get the images from suject (1, 2, 4)
folder = "subject" + num2str(subject); % Name of the folder
facial_expression = 1; % Select facial expression (1, 2, 3, 4, 5)

im_set_left = imageSet(convertStringsToChars(folder + "/" + folder + "_Left"));
im_set_middle = imageSet(convertStringsToChars(folder + "/" + folder + "_Middle"));
im_set_right = imageSet(convertStringsToChars(folder + "/" + folder + "_Right"));

im_left = im2double(read(im_set_left, facial_expression)); % Left image of the subject
im_middle = im2double(read(im_set_middle, facial_expression)); % Middle image of the subject
im_right = im2double(read(im_set_right, facial_expression)); % Right image of the subject

%% Camera calibration
% The images were taken with a stereo camera.
% For calibration purposes we have pairs of two pairs of images to get a
% more accurate description of the camera model. These are left-middle
% and middle-right.
% To get the intrisic and extrinsic parameters of the cameras we used the
% stereoCameraCalibrator app.

% Load in the parameters which were obtained from the app
load("Camera Data\calib2StereoParamsLM.mat");
load("Camera Data\calib2StereoParamsMR.mat");
load("Camera Data\calib2StereoParamsLR.mat");
% figure; showExtrinsics(stereoParamsLM); % show extrinsic parameters of LM
% figure; showExtrinsics(stereoParamsMR); % show extrinsic parameters of MR

%% Background remover
% Remove the background from the photo and extract only the face of the
% subject
face_left_mask = background_remover(im_left);
face_middle_mask = background_remover(im_middle);
face_right_mask = background_remover(im_right);

% Plot the three masked images (Left, middle, right)
figure;
subplot(1,3,1);imshow(face_left_mask.*im_left);title('Left');
subplot(1,3,2);imshow(face_middle_mask.*im_middle);title('Middle');
subplot(1,3,3);imshow(face_right_mask.*im_right);title('Right');
%% Stereo Rectification
% Reproject image planes onto a common plane parallel to the line between
% camera center
% Rectify the images of the subject in this order: left -> middle, 
% middle -> right
[im_left_rect, im_middle_left_rect] = rectifyStereoImages(...
    im_left, im_middle, calib2StereoParamsLM,'OutputView','full');
    
[im_middle_right_rect, im_right_rect] = rectifyStereoImages(...
    im_middle, im_right, calib2StereoParamsMR, 'OutputView','full');

% Do the same operation for the masks Left -> Middle -> Right
[mask_left_rect, mask_middle_left_rect] = rectifyStereoImages(...
    face_left_mask, face_middle_mask, calib2StereoParamsLM,'OutputView','full');
    
[mask_middle_right_rect, mask_right_rect] = rectifyStereoImages(...
    face_middle_mask, face_right_mask, calib2StereoParamsMR, 'OutputView','full');

[im_LR_l_rect, im_LR_r_rect, reprojection_mat_LR] = rectifyStereoImages(...
    im_left, im_right, calib2StereoParamsLR,'OutputView','full');

[im_LR_l_rect_mask, im_LR_r_rect_mask] = rectifyStereoImages(...
    face_left_mask, face_right_mask, calib2StereoParamsLR,'OutputView','full');

% Unify brightness
im_middle_left_rect = mimic_colorspace(im_middle_left_rect ,im_left_rect);
im_right_rect = mimic_colorspace(im_right_rect ,im_middle_right_rect);
im_LR_r_rect = mimic_colorspace(im_LR_r_rect ,im_LR_l_rect);

% Create rectified images cell arrays (for ease of access afterwards): 
% column1: first rectified image, column2: second rectified image
LM_rec = {im_left_rect im_middle_left_rect};   % Left-Middle
MR_rec = {im_middle_right_rect, im_right_rect}; % Middle-Right
LR_rec = {im_LR_l_rect, im_LR_r_rect};          % Left-Right

% Create rectified masks cell arrays (for ease of access afterwards): 
% column1: mask of first image, column2: mask of second image
LM_mask = {mask_left_rect, mask_middle_left_rect};      % Left-Middle
MR_mask = {mask_middle_right_rect, mask_right_rect};    % Middle-Right
LR_mask = {im_LR_l_rect_mask, im_LR_r_rect_mask};       % Left-Right

% Create rectified masked imgs cell arrays (for ease of access afterwards): 
% column1: first masked image, column2: second masked image
LM_masked = {LM_rec{1,1}.*LM_mask{1,1}, LM_rec{1,2}.*LM_mask{1,2}};% Left-Middle
MR_masked = {MR_rec{1,1}.*MR_mask{1,1}, MR_rec{1,2}.*MR_mask{1,2}};% Middle-Right
LR_masked = {LR_rec{1,1}.*LR_mask{1,1}, LR_rec{1,2}.*LR_mask{1,2}};% Left-Right

%% Stereo Matching and disparity maps
ref = 1; % dmap to be calculated wrt 1:first image, 2:second image 

% Minimum disparity for each subject was manually calculated. This was 
% done by inputing the Anaglyph to imtool, then measuring the distance
% between the deepest two corresponding points (of the neck for example).
% mdisp -> minimum disparity
mdisp_LM = [250 280 0 310]; %[subject1 subject2 dummy subject4]
mdisp_MR = [260 265 0 325]; %[subject1 subject2 dummy subject4]
mdisp_LR = [530 560 0 691]; %[subject1 subject2 dummy subject4]

% Disparity range can be then calculated. (maxD-minD) should be <= 128 and
% divisible by 8. So, to achieve max range, maxD = minD +128.
drange_LM = [mdisp_LM(subject), mdisp_LM(subject) + 128];
drange_MR = [mdisp_MR(subject), mdisp_MR(subject) + 128];
drange_LR = [mdisp_LR(subject), mdisp_LR(subject) + 128];

% Left-middle disparity map
[dmap_LM, reliability_map_LM]= calculateDisparityMap(LM_masked, drange_LM, ref);
% Left-middle enhance disparity map
[dmap_LM, reliability_map_LM] = enhanceDmap(dmap_LM, reliability_map_LM, ...
    LM_masked, drange_LM, LM_mask, ref);
% show_pc(dmap_LM, drange_LM,calib2StereoParamsLM, LM_masked{ref});

% Middle-right disparity map
[dmap_MR, reliability_map_MR]= calculateDisparityMap(MR_masked, drange_MR, ref);
% Middle-right enhance disparity map
[dmap_MR, reliability_map_MR] = enhanceDmap(dmap_MR, reliability_map_MR, ...
    MR_masked, drange_MR, MR_mask, ref);
% show_pc(dmap_MR, drange_MR,calib2StereoParamsMR, MR_masked{ref});

% Left-right disparity map
[dmap_LR, reliability_map_LR]= calculateDisparityMap(LR_masked, drange_LR, ref);
% Left-right enhance disparity map
[dmap_LR, reliability_map_LR] = enhanceDmap(dmap_LR, reliability_map_LR, ...
    LR_masked, drange_LR, LR_mask, ref);
% show_pc(dmap_LR, drange_LR,calib2StereoParamsLR, LR_masked{ref});

%% Create point clouds
[pc_LM, pc_LM_denoised] = Point_Clouds(dmap_LM, calib2StereoParamsLM, LM_masked{ref});
[pc_MR, pc_MR_denoised] = Point_Clouds(dmap_MR, calib2StereoParamsMR, MR_masked{ref});
[pc_LR, pc_LR_denoised] = Point_Clouds(dmap_LR, calib2StereoParamsLR, LR_masked{ref});
%% Point clouds merging
pc_LM_aligned = pc_LM_denoised; % For convention (aligned)

% % Align point clouds by using ICP algorithm
% [tform_MR_to_global,pc_MR_aligned,rmseMR_LM,originalrmseMR_LM] = pcregistericp_err(pc_MR_denoised,pc_LM_denoised,MaxIterations=150,verbose = true);
% [tform_LR_to_global,pc_LR_aligned,rmseLR_LM,originalrmseLR_LM] = pcregistericp_err(pc_LR_denoised,pc_LM_denoised,MaxIterations=150,verbose = true);

% Load saved point clouds (for time optimization)
pc = "PointClouds" + num2str(subject);
pc_MR_aligned = pcread(pc +"/pc_MR_aligned.ply");
pc_LR_aligned = pcread(pc +"/pc_LR_aligned.ply");

% Load save transformation matrices (for time optimization)
tforms = "tforms" + num2str(subject); % Name of the folder
load(tforms+"/tform_LR_to_global.mat");
load(tforms+"/tform_MR_to_global.mat");

% Show Point Cloud LM and new Point Cloud MR (aligned with Point Cloud LM)
figure; subplot(1,2,1);pcshow(pc_MR_aligned, BackgroundColor=[1,1,1]);
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
subplot(1,2,2);pcshow(pc_LM_aligned,BackgroundColor=[1,1,1]);
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
sgtitle("$$New$$ $$Point$$ $$Cloud$$ $$MR$$ $$and$$ $$Point$$ $$Cloud$$ $$LM$$",'interpreter','latex');
% Show Point Cloud LM and new Point Cloud LR (alligned with Point Cloud LR)
figure; subplot(1,2,1);pcshow(pc_LR_aligned, BackgroundColor=[1,1,1]);
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
subplot(1,2,2);pcshow(pc_LM_aligned, BackgroundColor=[1,1,1]);
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
sgtitle("$$New$$ $$Point$$ $$Cloud$$ $$LR$$ $$and$$ $$Point$$ $$Cloud$$ $$LM$$",'interpreter','latex');

%% Merge point clouds
merge_size = 0.02; % unit: meter

% Merging LM and MR
pc_global = pcmerge(pc_LM_aligned, pc_MR_aligned, merge_size);
figure;pcshow(pc_global, BackgroundColor=[1,1,1]);
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
title('Merged point clouds LM and MR', 'Interpreter','latex');
%% Create 3D surface meshes
[Tri_LM, img_LM_line] = make_3D_surface_mesh(dmap_LM, pc_LM, ...
    ~reliability_map_LM, LM_rec{ref});
title('3D Surface Mesh LM')

[Tri_MR, img_MR_line] = make_3D_surface_mesh(dmap_MR, pc_MR, ...
    ~reliability_map_MR, MR_rec{ref});
title('3D Surface Mesh MR')

[Tri_LR, img_LR_line] = make_3D_surface_mesh(dmap_LR, pc_LR, ...
    ~reliability_map_LR, LR_rec{ref});
title('3D Surface Mesh LR')
%% Merge 3D surface meshes
% Align each mesh to the left-middle 3D mesh
P_MR_aligned = transformPointsForward(tform_MR_to_global,Tri_MR.Points);
Tri_MR_aligned = triangulation(Tri_MR.ConnectivityList, P_MR_aligned);

P_LR_aligned = transformPointsForward(tform_LR_to_global,Tri_LR.Points);
Tri_LR_aligned = triangulation(Tri_LR.ConnectivityList, P_LR_aligned);

% Merge MR mesh to LM mesh to obtain LMMR mesh
merged_points = [Tri_LM.Points;Tri_MR_aligned.Points];
merged_conn = [Tri_LM.ConnectivityList;...
    Tri_MR_aligned.ConnectivityList+ size(Tri_LM.Points,1)];
Tri_LMMR = triangulation(merged_conn, merged_points);
merg_color_line_LMMR = [img_LM_line;img_MR_line];

% Merge LR mesh to LM mesh to obtain LMLR mesh
merged_points = [Tri_LM.Points;Tri_LR_aligned.Points];
merged_conn = [Tri_LM.ConnectivityList;...
    Tri_LR_aligned.ConnectivityList+ size(Tri_LM.Points,1)];
Tri_LMLR = triangulation(merged_conn, merged_points);
merg_color_line_LMLR = [img_LM_line;img_LR_line];

% Merge LR mesh to MR mesh to obtain MRLR mesh
merged_points = [Tri_MR_aligned.Points;Tri_LR_aligned.Points];
merged_conn = [Tri_MR_aligned.ConnectivityList;...
    Tri_LR_aligned.ConnectivityList+ size(Tri_MR_aligned.Points,1)];
Tri_MRLR = triangulation(merged_conn, merged_points);
merg_color_line_MRLR = [img_MR_line;img_LR_line];

% Merge LR mesh to LMMR mesh to obtain LMMRLR mesh
merged_points = [Tri_LMMR.Points;Tri_LR_aligned.Points];
merged_conn = [Tri_LMMR.ConnectivityList;...
    Tri_LR_aligned.ConnectivityList+ size(Tri_LMMR.Points,1)];
Tri_LMMRLR = triangulation(merged_conn, merged_points);
merg_color_line_LMMRLR = [merg_color_line_LMMR;img_LR_line];

% Visualize merged mesh
figure;
subplot(2,2,1);
TM = trimesh(Tri_LMMR);
title('Merged LMMR mesh');
set(TM,'FaceVertexCData',merg_color_line_LMMR); % set colors to input image
set(TM,'Facecolor','interp');
% set(TM,'FaceColor','red');            % if you want a colored surface
set(TM,'EdgeColor','none');             % suppress the edges
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
axis([-250 250 -250 250 400 900])
set(gca,'xdir','reverse')
set(gca,'zdir','reverse')
daspect([1,1,1])
axis tight

subplot(2,2,2);
TM = trimesh(Tri_LMLR);
title('Merged LMLR mesh');
set(TM,'FaceVertexCData',merg_color_line_LMLR); % set colors to input image
set(TM,'Facecolor','interp');
% set(TM,'FaceColor','red');            % if you want a colored surface
set(TM,'EdgeColor','none');             % suppress the edges
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
axis([-250 250 -250 250 400 900])
set(gca,'xdir','reverse')
set(gca,'zdir','reverse')
daspect([1,1,1])
axis tight

subplot(2,2,3);
TM = trimesh(Tri_MRLR);
title('Merged MRLR mesh');
set(TM,'FaceVertexCData',merg_color_line_MRLR); % set colors to input image
set(TM,'Facecolor','interp');
% set(TM,'FaceColor','red');            % if you want a colored surface
set(TM,'EdgeColor','none');             % suppress the edges
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
axis([-250 250 -250 250 400 900])
set(gca,'xdir','reverse')
set(gca,'zdir','reverse')
daspect([1,1,1])
axis tight

subplot(2,2,4);
TM = trimesh(Tri_LMMRLR);
title('Merged LMMRLR mesh');
set(TM,'FaceVertexCData',merg_color_line_LMMRLR); % set colors to input image
set(TM,'Facecolor','interp');
% set(TM,'FaceColor','red');            % if you want a colored surface
set(TM,'EdgeColor','none');             % suppress the edges
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
axis([-250 250 -250 250 400 900])
set(gca,'xdir','reverse')
set(gca,'zdir','reverse')
daspect([1,1,1])
axis tight


% %% Accuracy
% % Initial RMSE
% disp('RMSE before alignment')
% disp(originalrmseMR_LM(1))
% disp(originalrmseLR_LM(1))
% 
% 
% % After alignment
% disp('RMSE after alignment')
% disp(rmseMR_LM)
% disp(rmseLR_LM)

%% Robustness
% Load eye location in each point cloud
%[LM MR LR]
load("Robustness/subj1_eye_after.mat");
load("Robustness/subj1_eye_before.mat");
load("Robustness/subj2_eye_after.mat");
load("Robustness/subj2_eye_before.mat");
load("Robustness/subj4_eye_after.mat");
load("Robustness/subj4_eye_before.mat");
% Distances before
%[LM-MR MR-LR LM-LR]
subj1_distance_before = [...
    sqrt(sum((subj1_eye_before{1}-subj1_eye_before{2}).^2));
    sqrt(sum((subj1_eye_before{2}-subj1_eye_before{3}).^2));
    sqrt(sum((subj1_eye_before{1}-subj1_eye_before{3}).^2));];

subj2_distance_before = [...
    sqrt(sum((subj2_eye_before{1}-subj2_eye_before{2}).^2));
    sqrt(sum((subj2_eye_before{2}-subj2_eye_before{3}).^2));
    sqrt(sum((subj2_eye_before{1}-subj2_eye_before{3}).^2));];

subj4_distance_before = [...
    sqrt(sum((subj4_eye_before{1}-subj4_eye_before{2}).^2));
    sqrt(sum((subj4_eye_before{2}-subj4_eye_before{3}).^2));
    sqrt(sum((subj4_eye_before{1}-subj4_eye_before{3}).^2));];

% Distances after
%[LM-MR MR-LR LM-LR]
subj1_distance_after = [...
    sqrt(sum((subj1_eye_after{1}-subj1_eye_after{2}).^2));
    sqrt(sum((subj1_eye_after{2}-subj1_eye_after{3}).^2));
    sqrt(sum((subj1_eye_after{1}-subj1_eye_after{3}).^2));];

subj2_distance_after = [...
    sqrt(sum((subj2_eye_after{1}-subj2_eye_after{2}).^2));
    sqrt(sum((subj2_eye_after{2}-subj2_eye_after{3}).^2));
    sqrt(sum((subj2_eye_after{1}-subj2_eye_after{3}).^2));];

subj4_distance_after = [...
    sqrt(sum((subj4_eye_after{1}-subj4_eye_after{2}).^2));
    sqrt(sum((subj4_eye_after{2}-subj4_eye_after{3}).^2));
    sqrt(sum((subj4_eye_after{1}-subj4_eye_after{3}).^2));];
