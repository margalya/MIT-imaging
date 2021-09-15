function [res] = testTempFromWings(obj, pic, x0, y0)

firstGuess = [obj.OD obj.sigmaX obj.sigmaY obj.C];
lower = [0.2*obj.OD 0.1*obj.sigmaX 0.1*obj.sigmaY -1];
upper = [2*obj.OD 10*obj.sigmaX 10*obj.sigmaY 1 ];
% binW = appData.options.avgWidth;
[h, w] = size(pic);
X = (1 : w) + x0 - 1;
Y = (1 : h) + y0 - 1;
% nsigma = 0 : 0.25 : 2; %number of cloud sigmas to omit using weights
nsigma = 1.7;

for i = 1 : length(nsigma)
    % prepare weights matrix
    weights = ones(size(pic));
    [Xmask,Ymask] = meshgrid( obj.ROILeft:obj.ROIRight, obj.ROITop:obj.ROIBottom );
    mask = ((Xmask - obj.x0)/(obj.sigmaX*nsigma(i))).^2 + ((Ymask-obj.y0)/(obj.sigmaY*nsigma(i))).^2 < 1;
%     [Xmask,Ymask] = meshgrid( X, Y );
%     mask = ((Xmask - obj.x0)/(obj.sigmaX*nsigma(i))).^2 + ((Ymask - obj.y0)/(obj.sigmaY*nsigma(i))).^2 < 1;
    weights( mask ) = 0;
%     figure;
%     imagesc(weights);
%     axis equal
    
    [Xout, Yout, binnedPic, weightsOut] = prepareSurfaceData(X, Y, pic, weights);
    s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', firstGuess, 'Lower', lower, 'Upper', upper, 'Weights', weightsOut);
    f = fittype('OD*exp(-(x-x0)^2/2/sigmaX^2 - (y-y0)^2/2/sigmaY^2) + C', 'coefficients', {'OD', 'sigmaX', 'sigmaY', 'C'},...
        'problem', {'x0', 'y0'}, 'independent', {'x','y'}, 'dependent', 'z', 'options', s);
    res = fit([Xout, Yout], binnedPic, f, 'problem', {obj.res.x0, obj.res.y0});
    OD(i) = res.OD;
    sigmaX(i) = res.sigmaX; %#ok<*AGROW>
    sigmaY(i) = res.sigmaY; 
    C(i) = res.C;
    
    conf = confint(res);
    conf = (conf(2,:) - conf(1,:))/2;
    ODErr(i) = conf(1);
    sigmaXErr(i) = conf(2);
    sigmaYErr(i) = conf(3);
    CErr(i) = conf(4);

end

% figure;
% errorbar(nsigma, sigmaX, sigmaXErr, 'o')
% hold on;
% errorbar(nsigma, sigmaY, sigmaYErr, 'o')
% xlabel('Number of excluded cloud sigmas from center')
% ylabel('Size X, Size Y [pixels]')
% legend('Size X', 'Size Y')
% figure;
% errorbar(nsigma, C, CErr, 'o')
% xlabel('Number of excluded cloud sigmas from center')
% ylabel('C coefficient')
% figure;
% errorbar(nsigma, OD, ODErr, 'o')
% xlabel('Number of excluded cloud sigmas from center')
% ylabel('OD')

% Save result to file, concat
% fileName = 'E:\Dropbox (Personal)\MATLAB\MIT imaging New\sizeFromWings.m';
% m = matfile(fileName, 'Writable', true); %Note: writable is true by default IF the file does not exist
% if ~exist(fileName, 'file')
%     m.sigmaX = sigmaX;
%     m.sigmaY = sigmaY;
% else
%     m.sigmaX = [m.sigmaX sigmaX];
%     m.sigmaY = [m.sigmaY sigmaY];
% end
% clear m;
end
