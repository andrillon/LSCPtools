function smooth_vec = Gaussian_Smooth_TA2(original_vec, my_N)
    
if my_N~=0
    g = gausswin(1000,my_N);
    g = g ./ sum(g);
    
    %% Since this code might be used for EEG, we must rectify first !!
    original_vec = abs(original_vec);
%     taped_vector=[original_vec(floor(5/my_N):-1:1) original_vec original_vec(end:-1:end-floor(5/my_N)+1)];
    smooth_vec = conv(original_vec, g);
    smooth_vec = smooth_vec(floor(length(g)/2):end - floor(length(g)/2));
else
    smooth_vec=original_vec;
end