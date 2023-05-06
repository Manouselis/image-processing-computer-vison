function show_pc(Dmap, range,stereoParams,pointsColor)
    figure;
    imshow(Dmap,range);
    colormap jet;
    colorbar;

    xyzPoints = reconstructScene(Dmap,stereoParams);
    pc = pointCloud(xyzPoints,"Color",pointsColor);
    pc = pcdenoise(pc);
    % pc_LM = pcdownsample(pc_LM, 'nonuniformGridSample', 15);
    figure;
    pcshow(pc);
end