function [output] = LowPassFilter(input, filterSize, filterSigma )

%1/2D low pass filter for images, using a convolution with a Gaussian
%1 or 2 dimensions filtering is set by input size (vector / matrix)
if filterSigma<1
    filterSigma = 1;
end
if min(size(input))>1 %input is matrix
    filt = GaussianFilter2D( filterSize, filterSigma);
    output = conv2(double(input),filt,'same');
    nanMask = nan*ones(size(output)); %nanmask changes the edges of the pictures to nan, in order to avoid edge effects of the filter
    nanMask(filterSigma : end-filterSigma , filterSigma : end-filterSigma) = 1;
    output = nanMask.*output;
%     figure;imagesc(output)
elseif min(size(input))==1 %input is vector
    filt = GaussianFilter1D( filterSize, filterSigma);
    output = conv(double(input),filt,'same');
%     nanMask = nan*ones(size(output));
%     nanMask(filterSigma : end-filterSigma) = 1;
%     output = nanMask.*output;
    %     figure;plot(output)
else
    output = input;
    
end
end

function [h] = GaussianFilter2D(p2, std)
siz   = ([p2 p2]-1)/2;

[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
arg   = -(x.*x + y.*y)/(2*std*std);

h     = exp(arg);
h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));
if sumh ~= 0,
    h  = h/sumh;
end
end

function [h] = GaussianFilter1D(p2, std)
siz   = (p2-1)/2;

[x] = -siz(1):siz(1);
arg   = -(x.*x)/(2*std*std);

h     = exp(arg);
h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));
if sumh ~= 0,
    h  = h/sumh;
end
end
