normaly = @(mn, variance, x) 2.718283^(- (x - mn)*(x - mn)/(2*variance))/realsqrt(2*pi*variance);

sd = input('What is the sample standard error? ');
sd2 = sd*sd;
obtained = input('What is the sample mean? ');

uniform = input('is the distribution of p(population value|theory) uniform? 1= yes 0=no ');

if uniform == 0
    meanoftheory = input('What is the mean of p(population value|theory)? ');
    sdtheory = input('What is the standard deviation of p(population value|theory)? ');
    omega = sdtheory*sdtheory;
    tail = input('is the distribution one-tailed or two-tailed? (1/2) ');
end


if uniform == 1
    lower = input('What is the lower bound? ');
    upper = input('What is the upper bound? ');
end



area = 0;
if uniform == 1
    theta = lower;
else theta = meanoftheory - 5*(omega)^0.5;
end
if uniform == 1
    incr = (upper- lower)/2000;
else incr =  (omega)^0.5/200;
end

for A = -1000:1000
    theta = theta + incr;
    if uniform == 1
        dist_theta = 0;
        if and(theta >= lower, theta <= upper)
            dist_theta = 1/(upper-lower);
        end
    else %distribution is normal
        if tail == 2
            dist_theta = normaly(meanoftheory, omega, theta);
        else
            dist_theta = 0;
            if theta > 0
                dist_theta = 2*normaly(meanoftheory, omega, theta);
            end
        end
    end
    
    height = dist_theta * normaly(theta, sd2, obtained); %p(population value=theta|theory)*p(data|theta)
    area = area + height*incr; %integrating the above over theta
end


Likelihoodtheory = area
Likelihoodnull = normaly(0, sd2, obtained)
Bayesfactor = Likelihoodtheory/Likelihoodnull