function smooth_vec = Gaussian_Smooth_TA(original_vec, my_sigma)
    
if my_sigma~=0
    g = get_gauss(my_sigma);
    g = g ./ sum(g);
    
    %% Since this code might be used for EEG, we must rectify first !!
%     original_vec = abs(original_vec);
    taped_vector=[original_vec(floor(5*my_sigma):-1:1) original_vec original_vec(end:-1:end-floor(5*my_sigma)+1)];
    smooth_vec = conv(taped_vector, g);
    smooth_vec = smooth_vec(floor(length(g)/2)+floor(5*my_sigma):end - floor(length(g)/2)-floor(5*my_sigma));
else
    smooth_vec=original_vec;
end