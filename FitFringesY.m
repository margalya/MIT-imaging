classdef FitFringesY < FitTypes
%FITFRINGESY Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'FitFringesY';
   end
   properties (SetAccess = private )
       res = [];
       gof = [];
       output = [];
       a = -1;
       x0 = -1;
       w = -1;
       c = -1;
       conf = [];
   end

   methods
       function appData = analyze(obj, appData) % do the analysis
           % 1D fit
           fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussian};
           if ( isempty( fitObj.yRes ) )
               tmpFitType = appData.data.fitType;
               appData.data.fitType = appData.consts.fitTypes.oneDGaussian;
               appData = appData.data.fits{appData.consts.fitTypes.oneDGaussian}.analyze(appData);
               fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussian};
               appData.data.fitType = tmpFitType;
           end
           
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
           [h, w] = size(pic);
           
           % fit
           [xData, yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(...
               fitObj.xCenter, fitObj.yCenter, appData.options.avgWidth);
           xData = [1:length(yData)];
           yDataTemp=yData-(fitObj.yRes.Cy+fitObj.yRes.Ay*exp(-(xData-(fitObj.yRes.y0-y0)).^2./2./fitObj.yRes.sigmaY^2)); %substract gaussain from fringes, to make fft easier
           
           lambda = findWavelength( yDataTemp, fitObj.yRes.sigmaY );
           firstGuess = [fitObj.yRes.Ay, 0, 1.5*(max(yDataTemp)-min(yDataTemp)), fitObj.yRes.sigmaY, fitObj.yRes.y0-y0, lambda];
           lower([1 4 5 6]) = 0.7*firstGuess([1 4 5 6]); %a, w, x0
           upper([1 4 5 6]) = 1.3*firstGuess([1 4 5 6]); %a, w, x0
           lower(2) = -Inf; %phi
           upper(2) = Inf; %phi
           lower(3) = 0; %v
           upper(3) = 1; %v
           lower(6) = firstGuess(6)*0.9; %lambda
           upper(6) = firstGuess(6)*1.1; %lambda
%                       lower(6) = 25/2.3; %lambda
%                       upper(6) = 40/2.3; %lambda
           warning('off', 'MATLAB:class:cannotUpdateClass:Changed')
           
           fo = fitoptions('method','NonlinearLeastSquares','Robust','On', ...
               'DiffMinChange',1e-010,'DiffMaxChange',1e-005,'MaxFunEvals',4000, ...
               'MaxIter',8000,'TolFun',1e-010,'TolX',1e-010, 'Startpoint', firstGuess, ...
               'Lower', lower, 'Upper', upper);
           ft = fittype('a*exp(-(x-x0)^2/2/w^2)*(1+v*sin(2*pi/lambda*(x-length(x)/2)+phi))+c',...
                'dependent',{'y'},'independent',{'x'}, 'coefficients',{'a', 'phi', 'v', 'w', 'x0', 'lambda'}, 'problem', {'c'});
            % c is a problem parameter, in order to reduce error in determining te phase.
            
            % calcualte c from mean background level
%             tempPic = pic; tempPic( 1 : min(size(pic,1), round(fitObj.yRes.y0 - fitObj.yStart + 4.5*fitObj.yRes.sigmaY)),...
%                 round(fitObj.xRes.x0 - fitObj.xStart - 7.5*fitObj.xRes.sigmaX) : round(fitObj.xRes.x0 - fitObj.xStart + 7.5*fitObj.xRes.sigmaX) ) = nan;
% %             figure; imagesc(tempPic)
%             background = nanmean(tempPic(:));
%             clear tempPic
%             visibilityFFT( yData - fitObj.yRes.Cy )
%             ft = fittype([ num2str(obj.a) '*exp(-(x-' num2str(obj.x0) ')^2/2/' num2str(obj.w) '^2)*(1+v*sin(2*pi/lambda*(x-' num2str(obj.x0) ... 
%                 ')+phi))+' num2str(obj.c) ''], 'dependent',{'y'},'independent',{'x'}, 'coefficients',{ 'lambda', 'phi', 'v'});
           [obj.res, obj.gof, obj.output] = fit(xData',yData',ft,fo,'problem', {fitObj.yRes.Cy});
%            customFit4( obj, xData, yData )
%            customFit5( obj, xData, yData, appData.analyze.readDir, appData.save.picNo )
%            [visFFT,lambdaFFT]=visibility(xData', yData', 'Fourier');
%            assignin('base','visFFT',visFFT); assignin('base','lambdaFFT',lambdaFFT);
%            evalin('base','vis = [vis visFFT]; lambda = [lambda lambdaFFT];');
%            display(['Yoni''s calculation result: V = ' num2str(visFFT) '; lambda = ' num2str(lambdaFFT*appData.consts.cameras{appData.options.cameraType}.yPixSz*1e6) ' um;'])
           conf = confint(obj.res);
           obj.conf = (conf(2,:)-conf(1,:))/2;
%            obj.conf(6) = dlambda; %since lambda is not found by the fit, set its 'confidence bound' manually
           
            obj.a = obj.res.a;
            obj.x0 = obj.res.x0;
            obj.w = obj.res.w; %result is taken from fitFringesY, and not from 1DGaussian
            obj.c = obj.res.c;
            
           pxSz = appData.consts.cameras{appData.options.cameraType}.yPixSz;
           
            %set fit params
            obj.xCenter = round(fitObj.xRes.x0); % should be indexes (integers)
            obj.yCenter= round(fitObj.yRes.y0);
            obj.xUnitSize = fitObj.xRes.sigmaX;
            obj.yUnitSize = fitObj.yRes.sigmaY;
%             obj.maxVal = fitObj.yRes.Ay;
            obj.maxVal = obj.res.a *(1 + obj.res.v);

           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData);
           
            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
            [xData yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors( ...
                obj.xCenter, obj.yCenter, appData.options.avgWidth);

            obj.xData = xData;
            obj.xStart = x0;
            obj.yData = yData;
            obj.yStart = y0;

            % last
            appData.data.fits{appData.consts.fitTypes.fringesY} = obj;
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x);
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           normalizedROI = pic(y, x);
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-obj.c)/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = obj.c*ones(size(x));
           y = y-obj.yStart+1;
%            y = [y(1):0.01:y(end)];
           
%            yFit =  obj.res.a*exp(-(y-obj.res.x0).^2/2/obj.res.w^2) .*  ...
%                (1+obj.res.v*sin(2*pi*(y-obj.res.x0)/obj.res.lambda+obj.res.phi))+obj.res.c; %relative phase   
           yFit =  obj.res.a*exp(-(y-obj.res.x0).^2/2/obj.res.w^2) .*  ...
               (1+obj.res.v*sin(2*pi*(y-length(y)/2)/obj.res.lambda+obj.res.phi))+obj.res.c; %absolute phase
%            yFit =  obj.a*exp(-(y-obj.x0).^2/2/obj.w^2) .*  ...
%                (1+obj.res.v*sin(2*pi*(y-obj.x0)/obj.res.lambda+obj.res.phi))+obj.c;
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           cla(appData.ui.pFitResults)
           text(0, 1, ['\lambda = ' num2str(obj.res.lambda*appData.consts.cameras{appData.options.cameraType}.yPixSz*1e6)  ...
               ' +/- ' num2str(obj.conf(6)*appData.consts.cameras{appData.options.cameraType}.yPixSz*1e6) ' [\mum]'], 'Parent', appData.ui.pFitResults);
           text(0, 0.9, ['\phi = ' num2str(mod(obj.res.phi, 2*pi)) ' +/- ' num2str(obj.conf(2) ) ' [rad]' ], 'Parent', appData.ui.pFitResults);
           text(0, 0.8, ['visibility = ' num2str(obj.res.v) ' +/- ' num2str(obj.conf(3) )], 'Parent', appData.ui.pFitResults);
%            text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
%            text( 10, 160, 'fit function = A * {\ite}^{-(x-x_0)^2 / 2\sigma_x^2 - (y-y_0)^2 / 2\sigma_y^2} + C', 'fontsize', 12);
%            text( 50, 135, ['A = ' num2str(obj.OD)], 'fontsize', 12);
%            text( 50, 100, {['x_0 = ' num2str(obj.x0 * appData.data.camera.xPixSz * 1000) ' mm'], ...
%                ['\sigma_x = ' num2str(obj.sigmaX * appData.data.camera.xPixSz * 1000) ' mm']}, ...
%                'fontsize', 12);
%            text( 200, 100, {['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...
%                ['\sigma_y = ' num2str(obj.sigmaY * appData.data.camera.yPixSz * 1000) ' mm']}, ...
%                'fontsize', 12);
           text( 0, 0.7, ['c = ' num2str(obj.c)], 'Parent', appData.ui.pFitResults);
       end    
      
   end
end 
      
function [ lambda, phi] = findWavelength( yData, sigmaY )
% Finds the interference wavelength, in pixels, using fft
% Based on matlab's example for fft command

T = 1; %Sampling distance is 1 pixel
Fs = 1/T; %Sampling frequency
L = length(yData); %Length of signal
NFFT = 50*2^nextpow2(L); %50 factor increases the precision of the fft
Y = fft(yData,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1); %frequency vector
% figure;plot(f,2*abs(Y(1:NFFT/2+1)))
Yeff = 2*abs(Y(1:NFFT/2+1));
Yeff(f<(1/sigmaY/2)) = 0; % cut low frequencies - to cut the gaussian contribution to the fft signal
% figure;plot(f,Yeff) % plot after frequency cutoff
% xlabel('Frequency (Hz)')
[~,I] = max(Yeff); %find maximum frequency, use low freuquency cutoff, cutting spatial frequencies which are smaller than 1/sigmaY/2, meaning there are at least two fringes on gaussian
lambda=1/f(I);
% dlambda = (1/f(I-1) - 1/f(I+1))/2; %calculate error in lambda, using the two closest points to the maximum.
[~,Iraw] = max(abs(Y(1:end/2)));
phi = angle(Y(Iraw));
end