function my_gauss = get_gauss(my_sigma);

 
    if (nargin < 1), my_sigma = 300;, end;
 
    time_window = 8 * my_sigma;
    mu = round(time_window/2); 
    
     t = 1:time_window;

    my_gauss = exp(-((t-mu).^2)/(2*my_sigma*my_sigma)) ;