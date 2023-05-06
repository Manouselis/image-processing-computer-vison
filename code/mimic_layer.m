function adjustedlayer = mimic_layer(img,reference_img)
mean_reference = mean2(reference_img);
std_reference = std2(reference_img);
mean_img= mean2(img);
std_img= std2(img);

sigma_desired = std_reference;
mean_desired = mean_reference;

a = sigma_desired/std_img;
b = mean_desired - a*mean_img;

adjustedlayer = a*img + b;

end