
function binnedData = binning ( pic, binningFactor)
if binningFactor == 0
    binningFactor = 1;
end
[height width] = size(pic);
binnedData = zeros(floor(height / binningFactor), floor(width / binningFactor));

for i=1 : binningFactor
    for j=1 : binningFactor
        binnedData = binnedData + pic(i : binningFactor : floor(height / binningFactor)*binningFactor+ i - binningFactor, ...
            j : binningFactor : floor(width / binningFactor)*binningFactor + j - binningFactor);
    end
end
binnedData = binnedData ./ binningFactor^2;
end
