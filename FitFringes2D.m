classdef FitFringes2D < FitTypes
    %FITFRINGESY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant = true )
        ID = 'FitFringes2D';
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
            fitObj = appData.data.fits{appData.consts.fitTypes.fringesY};
            if ( isempty( fitObj.res ) )
                tmpFitType = appData.data.fitType;
                appData.data.fitType = appData.consts.fitTypes.oneDGaussian;
                appData = appData.data.fits{appData.consts.fitTypes.fringesY}.analyze(appData);
                fitObj = appData.data.fits{appData.consts.fitTypes.fringesY};
                appData.data.fitType = tmpFitType;
            end
            
            [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
            [h, w] = size(pic);
            
            % {'a', 'x0', 'y0', 'sigmaX', 'sigmaY', 'v', 'lambda', 'phi', 'c'}
            firstGuess = [fitObj.res.a appData.data.fits{appData.consts.fitTypes.oneDGaussian}.xRes.x0 fitObj.res.x0 + y0...
                appData.data.fits{appData.consts.fitTypes.oneDGaussian}.xRes.sigmaX fitObj.res.w fitObj.res.v fitObj.res.lambda fitObj.res.phi fitObj.res.c];
            lower = 0.95*firstGuess;
            upper = 1.05*firstGuess;
            lower(8) = -Inf; %phi
            upper(8) = Inf; %phi
            lower(7) = 0.98 * firstGuess(7); %lambda
            upper(7) = 1.02 * firstGuess(7); %lambda
            lower(9) = firstGuess(9) - abs(firstGuess(9)); %c
            upper(9) = firstGuess(9) + abs(firstGuess(9)); %c
            warning('off', 'MATLAB:class:cannotUpdateClass:Changed')
            
            fo = fitoptions('method','NonlinearLeastSquares','Robust','On', ...
                'DiffMinChange',1e-010,'DiffMaxChange',1e-005,'MaxFunEvals',1000, ...
                'MaxIter',4000,'TolFun',1e-010,'TolX',1e-010, 'Startpoint', firstGuess, ...
                'Lower', lower, 'Upper', upper);
            ft = fittype('a * exp(-(x-x0)^2/sigmaX^2/2) * (1+v*sin(2*pi/lambda*(y-length(y)/2)+phi)) * exp(-(y-y0)^2/sigmaY^2/2) + c',...
                'dependent',{'z'},'independent',{'x','y'}, 'coefficients',{'a', 'x0', 'y0', 'sigmaX', 'sigmaY', 'v', 'lambda', 'phi', 'c'});
            
            % 2D fit
            binW = appData.options.avgWidth;
%             if ( numel(pic) > 350^2 )
%                 binnedPic = binning(pic, binW);
%                 [h w] = size(binnedPic);
%                 [X, Y] = meshgrid([1 : binW : binW*w] +x0-1, [1 : binW : binW*h] +y0-1);
%                 [obj.res, obj.gof, obj.output] = fit( [X(:), Y(:)],pic(:), ft, fo);
%             else
                [X, Y] = meshgrid([1  : w] +x0-1, [1  : h] +y0-1);
                %                [obj.res, obj.gof, obj.output] = fit( [X(:), Y(:)],pic(:), ft, fo);
                [fitRes, fval, exitflag, output] = ...
                    fminsearch(@(p) fitFringes2D_scalar( p, X, Y, pic), firstGuess, optimset('TolX',1e-4, 'MaxIter', 10000, 'MaxFunEvals', 10000) );
%             end
            
            % visualize fit
            fitSurf = fitRes(1) * exp(-(X-fitRes(2)).^2./fitRes(4).^2./2) .* ( 1 + fitRes(6).*sin(2*pi/fitRes(7).*(Y-length(Y)./2) + fitRes(8) )) .* exp(-(Y-fitRes(3)).^2./fitRes(5).^2./2) +fitRes(9);
            figure; mesh(X,Y,fitSurf); hold on; stem3(X,Y,LowPassFilter(pic,30,1), 'marker', '.', 'LineStyle', 'None')
            figure; mesh(X,Y,fitSurf); hold on; stem3(X,Y,pic, 'marker', '.', 'LineStyle', 'None')
            %             figure; surf(X,Y,fitSurf); shading interp; hold on; scatter3(X(1:10:end),Y(1:10:end),pic(1:10:end),'marker','.');
            %             figure; surf(X,Y,pic); shading interp; hold on; mesh(X,Y,fitSurf);
            %             figure; surf(X,Y,fitSurf); hold on; mesh(X,Y,pic); shading interp;
            
            
            x = [1  : w] + x0 - 1;
            y = [1  : h] + y0 - 1;
            figure; plot(y, pic(:,round(fitRes(2) - x0))); hold on; plot(y, fitRes(1) * exp(-(0).^2./fitRes(4).^2./2) .* ( 1 + fitRes(6).*sin(2*pi/fitRes(7).*(y-length(y)./2) + fitRes(8) )) .* exp(-(y-fitRes(3)).^2./fitRes(5).^2./2) +fitRes(9), 'r');
            figure; plot(y, pic(:,round(fitRes(2) - fitRes(4) - x0))); hold on; plot(y, fitRes(1) * exp(-(- fitRes(4)).^2./fitRes(4).^2./2) .* ( 1 + fitRes(6).*sin(2*pi/fitRes(7).*(y-length(y)./2) + fitRes(8) )) .* exp(-(y-fitRes(3)).^2./fitRes(5).^2./2) +fitRes(9), 'r');
            figure; plot(y, pic(:,round(fitRes(2) + fitRes(4) - x0))); hold on; plot(y, fitRes(1) * exp(-(+ fitRes(4)).^2./fitRes(4).^2./2) .* ( 1 + fitRes(6).*sin(2*pi/fitRes(7).*(y-length(y)./2) + fitRes(8) )) .* exp(-(y-fitRes(3)).^2./fitRes(5).^2./2) +fitRes(9), 'r');
            figure; plot(x, pic(round(fitRes(3) - y0),:)); hold on; plot(x, fitRes(1) * exp(-(x-fitRes(2)).^2./fitRes(4).^2./2) .* ( 1 + fitRes(6).*sin(2*pi/fitRes(7).*(fitRes(3)-length(fitRes(3))./2) + fitRes(8) )) .* exp(-(0).^2./fitRes(5).^2./2) +fitRes(9), 'r');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 24/1/17 - end of work
            % until this point. The 2D fit more or less converges, but not sure about the quality
            % (visualization is hard).
            % I had to use 'fminsearch' instead of 'fit', it is easier.
            
            %            [obj.res, obj.gof, obj.output] = fit( [X(:), Y(:)],pic(:), ft, fo)
            conf = confint(obj.res);
            obj.conf = (conf(2,:)-conf(1,:))/2;
            
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
            
            text(10, 190, ['\lambda = ' num2str(obj.res.lambda*appData.consts.cameras{appData.options.cameraType}.yPixSz*1e6)  ...
                ' +/- ' num2str(obj.conf(6)*appData.consts.cameras{appData.options.cameraType}.yPixSz*1e6) ' [\mum]']);
            text(10, 170, ['\phi = ' num2str(mod(obj.res.phi, 2*pi)) ' +/- ' num2str(obj.conf(2) )]);
            text(10, 150, ['visibility = ' num2str(obj.res.v) ' +/- ' num2str(obj.conf(3) )]);
            %            text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
            %            text( 10, 160, 'fit function = A * {\ite}^{-(x-x_0)^2 / 2\sigma_x^2 - (y-y_0)^2 / 2\sigma_y^2} + C', 'fontsize', 12);
            %            text( 50, 135, ['A = ' num2str(obj.OD)], 'fontsize', 12);
            %            text( 50, 100, {['x_0 = ' num2str(obj.x0 * appData.data.camera.xPixSz * 1000) ' mm'], ...
            %                ['\sigma_x = ' num2str(obj.sigmaX * appData.data.camera.xPixSz * 1000) ' mm']}, ...
            %                'fontsize', 12);
            %            text( 200, 100, {['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...
            %                ['\sigma_y = ' num2str(obj.sigmaY * appData.data.camera.yPixSz * 1000) ' mm']}, ...
            %                'fontsize', 12);
            text( 10, 130, ['c = ' num2str(obj.c)]);
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

function ret = fitFringes2D_scalar( p, X, Y, pic )
% {'a', 'x0', 'y0', 'sigmaX', 'sigmaY', 'v', 'lambda', 'phi', 'c'}
% a = p(1);
% x0 = p(2);
% y0 = p(3);
% sigmaX = p(4);
% sigmaY = p(5);
% v = p(6)
% lambda = p(7)
% phi = p(8)
% c = p(9)

ret = p(1) * exp(-(X-p(2)).^2./p(4).^2./2) .* ( 1 + p(6).*sin(2*pi/p(7).*(Y-length(Y)./2) + p(8) )) .* exp(-(Y-p(3)).^2./p(5).^2./2) +p(9) - pic;
% ret = p(1)*( exp( -0.5 * (X - p(2)).^2./p(4).^2 - 0.5 * (Y - p(3)).^2./p(5).^2 ) ) + p(6) - pic;
ret = sum(sum(ret.^2));
end