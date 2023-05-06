function adjustedimage = mimic_colorspace(img,reference_img)
adjustedimage = zeros(size(img),'like',img);
for i=1:1:3
    adjustedimage(:,:,i) = mimic_layer(img(:,:,i),reference_img(:,:,i));
end
end