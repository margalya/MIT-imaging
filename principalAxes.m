function  principalAxes(appData)
% Find principle axes of a given image by rotation and fit

pic = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
theta = 0:10:180;
sigmaX = zeros(size(theta));
SigmaXErr = zeros(size(theta));
sigmaY = zeros(size(theta));
SigmaYErr = zeros(size(theta));
fig = figure;
for i = 1 : length(theta)
    fitPic  = pic + 10;
    fitPic = imrotate(fitPic, theta(i), 'bilinear', 'crop'); %'nearest'
    fitPic(fitPic==0) = NaN; %ignore zero values from the rotation
    fitPic = fitPic - 10;
    %     binnedData = binning ( pic, appData.options.avgWidth*2+1);
    binnedData = LowPassFilter(fitPic, 30, appData.options.avgWidth);
    [maxes, indexes] = max(binnedData); % find maximum
    [~, xPosMax] = max(maxes);
    yPosMax = indexes(xPosMax);
    xData = mean( fitPic(yPosMax-20 : yPosMax+20, :), 1);
    xData = xData(~isnan(xData));
    yData = mean( fitPic(:, xPosMax-20 : xPosMax+20), 2)';
    yData = yData(~isnan(yData));
    x = 1:length(xData);
    y = 1:length(yData);
    %     xData = LowPassFilter(xData, 15, appData.options.avgWidth);
    %     yData = LowPassFilter(yData, 15, appData.options.avgWidth);
    [ startPoint, lower, upper ] = calcInitialGuess( x, xData, appData.options.avgWidth );
    s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
    f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
    xRes = fit(x', xData', f);
    xBounds = calcBounds(xRes);
    [ startPoint, lower, upper ] = calcInitialGuess( y, yData, appData.options.avgWidth );
    s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
    f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
    yRes = fit(y', yData', f);
    yBounds = calcBounds(yRes);
    sigmaX(i) = xRes.sigmaX;
    SigmaXErr(i) = xBounds(3);
    sigmaY(i) = yRes.sigmaY;
    SigmaYErr(i) = yBounds(3);
    verticaShift = 0.1;
    figure(fig);
    subplot(2,1,1);
%     cla;
    plot(x, xData + verticaShift*i);
    hold on;
    plot(x, xRes(x) + verticaShift*i,'r');
    title('X axis fit')
    xlabel('Position [pixels]')
    ylabel('OD + shift')
    subplot(2,1,2)
%     cla;
    plot(y, yData + verticaShift*i);
    hold on;
    plot(y, yRes(y) + verticaShift*i,'r');
    title ('Y axis fit');
    xlabel('Position [pixels]')
    ylabel('OD + shift')
end

figure;
errorbar(theta, sigmaX, SigmaXErr, 'o')
hold on;
errorbar(theta, sigmaY, SigmaYErr, 'or')
legend({'Sigma X','Sigma Y'})
xlabel('Rotation angle [degrees]')
ylabel('Cloud size X,Y [pixels]')
figure;
plot(theta, sigmaX./sigmaY, 'o')
xlabel('Rotation angle [degrees]')
ylabel('SigmaX/SigmaY')
end

function  [ startPoint, lower, upper ] = calcInitialGuess( x, xData, filterSigma )
if filterSigma==0
    xData = xData( round(1.5): (end-round(1.5)) ); % omit points which are at the edges, to avoid filter edge effects
else
    xData = xData( round(1.5*filterSigma): (end-round(1.5*filterSigma)) ); % omit points which are at the edges, to avoid filter edge effects
end
% startPoint =[]; % {'Ax', 'x0', 'sigmaX', 'Cx'}
%                xDataS = smooth(xData, round(length(xData)/5) )'; %
%                smoothing dropped, since data is passed through a low pass
%                filter
[~, Ix] = max(xData);
CStartPoint = mean( xData(xData<(mean(xData)+std(xData))) ); %background start point is the mean value of all points which are no further than 1 sigma from mean value of xData
AStartPoint = max(xData) - CStartPoint; %Gaussian amplitude
xData(xData<(mean(xData)+std(xData))) = 0; % set points outside of peak to zero - to evaluate sigma
sigmaStartPoint = find(xData, 1, 'last' )-find(xData, 1 );

startPoint = [ AStartPoint x(Ix)+round(1.5*filterSigma) sigmaStartPoint CStartPoint]; % [Amplitude, x0, sigma, background]
lower(1) = 0.5 * startPoint(1);
upper(1) = 2 * startPoint(1);
lower(2) = startPoint(2) - 40;
upper(2) = startPoint(2) + 40;
lower(3) = 0.25 * startPoint(3);
upper(3) = 2*startPoint(3);
lower(4) = -2; % startPoint(4) - abs(startPoint(4));
upper(4) = 2; %startPoint(4) + abs(startPoint(4));

end

