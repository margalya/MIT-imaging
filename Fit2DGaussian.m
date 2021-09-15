classdef Fit2DGaussian < FitTypes
%FIT21DGAUSSIAN Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'Fit2DGaussian';
   end
   properties (SetAccess = private )
       res = [];
       gof = [];
%        output = [];
       conf = [];
       OD = -1;
       x0 = -1;
       y0 = -1;
       sigmaX = -1;
       sigmaY = -1;
       sigmaXWings = -1; %Wings only fit result of size, sigma>1.7
       sigmaYWings = -1;
       ODWings = -1;
       C = -1;
       calibMeanCountWithout = -1; % calibrated mean absortption photons count of without atoms image (in entire image, not ROI), calibrated to the withAtoms image level using C parameter from the fit
       photonCount = -1; %Fluorescence photon count, used for large signal images where the background is jittering.

   end

   methods
       function appData = analyze(obj, appData) % do the analysis
           % 1D fit
           fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussian};
           if ( isempty( fitObj.xRes ) )
               tmpFitType = appData.data.fitType;
               appData.data.fitType = appData.consts.fitTypes.oneDGaussian;
               appData = appData.data.fits{appData.consts.fitTypes.oneDGaussian}.analyze(appData);
               fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussian};
               appData.data.fitType = tmpFitType;
           end
           
%            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);  % the whole pic
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
           [h w] = size(pic);
           
            % 2D fit using matlab fit function
           firstGuess = [0.5*(fitObj.xRes.Ax+fitObj.yRes.Ay) fitObj.xRes.x0 fitObj.yRes.y0 ...
               fitObj.xRes.sigmaX fitObj.yRes.sigmaY 0.5*(fitObj.xRes.Cx+fitObj.yRes.Cy)];
           lower = [0 x0 y0 0 0 -30]; %#ok<PROPLC>
           upper = [10000 x0+w y0+h w h 30 ];%#ok<PROPLC>
           binW = appData.options.avgWidth;
           s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', firstGuess, 'Lower', lower, 'Upper', upper);
           f = fittype('OD*exp(-(x-x0)^2/2/sigmaX^2 - (y-y0)^2/2/sigmaY^2) + C', 'coefficients', {'OD', 'x0', 'y0', 'sigmaX', 'sigmaY', 'C'},...
               'independent', {'x','y'}, 'dependent', 'z', 'options', s);
            if ( numel(pic) > 450^2 )
                binnedPic = binning(pic, binW);
                [h, w] = size(binnedPic);
                X = (1 : binW : binW*w) +x0-1; %#ok<PROPLC>
                Y = (1 : binW : binW*h) +y0-1; %#ok<PROPLC>
                [X, Y, binnedPic] = prepareSurfaceData(X, Y, binnedPic);
                [obj.res, obj.gof, output] = fit([X, Y], binnedPic, f);
            else
                X = (1 : w) +x0-1; %#ok<PROPLC>
                Y = (1 : h) +y0-1; %#ok<PROPLC>
                [X, Y, binnedPic] = prepareSurfaceData(X, Y, pic);
                [obj.res, obj.gof, output] = fit([X, Y], binnedPic, f);
            end
            obj.conf = confint(obj.res);
            obj.conf = (obj.conf(2,:)-obj.conf(1,:))/2;
            
            % set 2D params
            obj.OD = obj.res.OD;
            obj.x0 = obj.res.x0;
            obj.y0 = obj.res.y0;
            obj.sigmaX = obj.res.sigmaX;
            obj.sigmaY = obj.res.sigmaY;
            obj.C = obj.res.C;
            
%             % set light level into saveParamVal
%             calibratedAtom = (appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic + 1)*exp(-obj.C)-1; % background picture, calibrated to the atoms level using C parameter from the fit
%             appData.save.saveParamVal = mean(calibratedAtom(:));
%             set(appData.ui.etParamVal, 'String', appData.save.saveParamVal)
            
            %%%%%%%%%%% calculate residuals of the 2D fit:
%             disp(['Mean 2D Gaussian fit resiuals RMS = ' num2str(sqrt(mean(output.residuals(:).^2))) ])
%             figure;
%              imagesc(LowPassFilter(reshape(output.residuals,size(pic)),30,5) )
% %             imagesc(reshape(output.residuals,size(pic)) )
%             colorbar
%             title('2D Gaussian fit residuals')
            %%%%%%%%%%%
            %set fit params
            obj.xCenter = round(obj.x0); % should be indexes (integers)
            obj.yCenter= round(obj.y0);
            obj.xUnitSize = obj.sigmaX;
            obj.yUnitSize = obj.sigmaY;
            obj.maxVal = obj.OD;
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData);
%            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           
%             %assign sizes and C coefficienct from wings fit - overwrites the full 2D result!
%             resWings = testTempFromWings(obj, pic, x0, y0);
%             obj.sigmaXWings = resWings.sigmaX;
%             obj.sigmaYWings = resWings.sigmaY;
%             obj.ODWings = resWings.OD;
%             %%% Overwrite regular results:
%             obj.sigmaX = resWings.sigmaX;
%             obj.sigmaY = resWings.sigmaY;
%             obj.xUnitSize = obj.sigmaX;
%             obj.yUnitSize = obj.sigmaY;
% % %             obj.C = resWings.C;
%             obj.OD = resWings.OD;
%             obj.maxVal = obj.OD;
            %%%
            if isfield(appData.consts.cameras{appData.options.cameraType}, 'photonPerADU') %If camera is setup for fluorescence imaging
                [obj.atomsNo, obj.photonCount] = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
                    [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
            else
                obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
                    [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
            end
%            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
%                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
            [xData yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(obj.xCenter, obj.yCenter, binW);
            
            obj.xData = xData;
            obj.xStart = x0;
            obj.yData = yData;
            obj.yStart = y0;
            
            obj.calibMeanCountWithout = (appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic + 1)*exp(-obj.C)-1; % background picture, calibrated to the withAtoms image level using C parameter from the fit, ROI only
%             obj.calibMeanCountWithout = obj.calibMeanCountWithout( obj.ROITop : obj.ROIBottom, obj.ROILeft : obj.ROIRight );
            obj.calibMeanCountWithout = mean2(obj.calibMeanCountWithout); %calibrated mean photon count in background picture, only in ROI
            % un-calibrated background photon count, for comparison:
%             test = mean2(appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic( obj.ROITop : obj.ROIBottom, obj.ROILeft : obj.ROIRight ));
%             disp(['Calibrated background photon count: ' num2str(obj.calibMeanCountWithout)])
%             disp(['Un-calibrated background photon count: ' num2str(test)])
            % last
            appData.data.fits{appData.consts.fitTypes.twoDGaussian} = obj;
%            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - obj.C;
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
           normalizedROI = obj.OD * exp( -0.5 * ((X-obj.x0).^2 / obj.sigmaX^2 + (Y-obj.y0).^2 / obj.sigmaY^2) );
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-obj.C)/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = obj.OD * exp( -0.5 * (x-obj.x0).^2 / obj.sigmaX^2 ) + obj.C;
           yFit = obj.OD * exp( -0.5 * (y-obj.y0).^2 / obj.sigmaY^2 ) + obj.C;
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           chipStart = appData.data.camera.chipStart;
           cla(appData.ui.pFitResults)
%            disp(obj.C)
           %            text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
           if isfield(appData.consts.cameras{appData.options.cameraType}, 'photonPerADU') %If camera is setup for fluorescence imaging, display number of photons
               %                photonCount = round(sum(appData.data.plots{2}.pic(:)) / 0.55 *0.28);
               text( 0, 1, ['Photon count: ' num2str(obj.photonCount)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           else %Display number of atoms
               text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           end
           text( 0, 0.8, 'fit function = A * {\ite}^{-(x-x_0)^2 / 2\sigma_x^2 - (y-y_0)^2 / 2\sigma_y^2} + C',...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.65,['A = ' num2str(obj.OD)], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.45, { ['x_0 = ' num2str(obj.x0 * appData.data.camera.xPixSz * 1000) ' mm;'], ...
               ['\sigma_x = ' num2str(obj.sigmaX * appData.data.camera.xPixSz * 1e6) ' um;']}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           
           text( 0.3, 0.45, {['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'],...
               ['\sigma_y = ' num2str(obj.sigmaY * appData.data.camera.yPixSz * 1e6) ' um']},...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           if isempty(obj.conf) %backwards compatability - before using fit in 2DGaussian
               text( 0, 0.25,['C = ' num2str(obj.C, '%.4f')], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
               AR = obj.sigmaX./obj.sigmaY;
               text( 0, 0.15,['X/Y size aspect ratio = ' num2str(AR)], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
           else
               text( 0, 0.25,['C = ' num2str(obj.C, '%.4f') ' +/- ' num2str(obj.conf(6), '%.4f')], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
               AR = obj.sigmaX./obj.sigmaY;
               ARErr = AR.*sqrt( (obj.conf(4)./obj.sigmaX).^2 + (obj.conf(5)./obj.sigmaY).^2 );
               text( 0, 0.15,['X/Y size aspect ratio = ' num2str(AR) ' +/- ' num2str(ARErr, '%.4f')], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
               meanWith = round( mean(appData.data.plots{appData.consts.plotTypes.withAtoms}.pic(:)) * 10) ./ 10 ;
               meanWithout = round( mean(appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic(:)) * 10) ./ 10;
               text( 0, 0.05,  ['Mean pixel count with, without = ' num2str(meanWith) ', ' num2str(meanWithout)], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
               tof = 3; %TOF in ms
               text( 0, -0.05,  ['Far field temperature (TOF = ' num2str(tof) ' ms) = ' num2str( round(singleShotTemperature(appData, obj, tof)*1e8)./1e2) ' uK'], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
           end
           end
   end
end

function [TMean] = singleShotTemperature(appData, obj, tof)

switch appData.options.atomType
    case 1
        mass = appData.consts.MNa;
    case 2
        mass = appData.consts.MLi6;
end
kB = appData.consts.Kb;

sx = obj.sigmaX  * appData.data.camera.xPixSz; % change to meters
sy = obj.sigmaY  * appData.data.camera.yPixSz; % change to meters

Tx = (sx./tof).^2 *mass / kB *1e6; %temperature in uK
Ty = (sy./tof).^2 *mass / kB *1e6;
TMean = (2*Tx + Ty)./3;

end