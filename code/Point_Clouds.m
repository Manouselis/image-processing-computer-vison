function [ptCloud_pair1, ptCloud_pair1_new] = Point_Clouds(disparityMap_pair1,reprojectionMatrix_pair1, im_rect)
point3d_pair1 = reconstructScene(disparityMap_pair1,reprojectionMatrix_pair1); % Based on whether we return reprojection Matrix from Stero Rectification (i.e., rectifyStereoImages)
%or reconstructScene(disparityMap_pair1,StereoParameters_pair1)
ptCloud_pair1 = pointCloud(point3d_pair1,"Color",im_rect);

% Remove noise (outliers) from pointCloud
%PreserveStructure=true in order to keep the correct dimension. Correct
%dimensions will be needed in a subsequent function provided by IPCV
%lectures that creates the 3D meshes
ptCloud_pair1_new = pcdenoise(ptCloud_pair1,PreserveStructure=true); 
figure; subplot(1,2,1);pcshow(ptCloud_pair1, BackgroundColor=[1,1,1]);% Compare before and after noise removal
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
subplot(1,2,2);pcshow(ptCloud_pair1_new, BackgroundColor=[1,1,1]);
xlabel('x (mm)','interpreter','latex');ylabel('y (mm)','interpreter','latex');zlabel('z (mm)','interpreter','latex');
sgtitle("$$Before$$ $$and$$ $$after$$ $$noise$$ $$removal$$",'interpreter','latex');
end



 