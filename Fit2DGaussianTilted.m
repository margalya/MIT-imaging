classdef Fit2DGaussianTilted < FitTypes
%FIT21DGAUSSIAN Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'Fit2DGaussianTilted';
   end
   properties (SetAccess = private )
       OD = -1;
       x0 = -1;
       y0 = -1;
       sigmaX = -1;
       sigmaY = -1;
       C = -1;
       
       fitRes = [];
       gof = [];
       output = [];
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
           [pic, x0, y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); %#ok<PROPLC> % only ROI
           [h, w] = size(pic);
           
           % 2D fit
           firstGuess = [0.5*(fitObj.xRes.Ax+fitObj.yRes.Ay) fitObj.xRes.x0 fitObj.yRes.y0 ...
               fitObj.xRes.sigmaX fitObj.yRes.sigmaY 0.5*(fitObj.xRes.Cx+fitObj.yRes.Cy)];
           lower = [0 x0 y0 0 0 -10]; %#ok<PROPLC>
           upper = [10 x0+w y0+h w h 10 ];%#ok<PROPLC>
           binW = appData.options.avgWidth;
           s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', firstGuess, 'Lower', lower, 'Upper', upper);
           f = fittype('OD*exp(-(x-x0)^2/2/sigmaX^2 - (y-y0)^2/2/sigmaY^2) + C', 'coefficients', {'OD', 'x0', 'y0', 'sigmaX', 'sigmaY', 'C'},...
               'independent', {'x','y'}, 'dependent', 'z', 'options', s);
            if ( numel(pic) > 350^2 )
                binnedPic = binning(pic, binW);
                [h, w] = size(binnedPic);
                X = (1 : binW : binW*w) +x0-1; %#ok<PROPLC>
                Y = (1 : binW : binW*h) +y0-1; %#ok<PROPLC>
                [X, Y, binnedPic] = prepareSurfaceData(X, Y, binnedPic);
                [obj.fitRes, obj.gof, obj.output] = fit([X, Y], binnedPic, f);
            else
                X = (1 : w) +x0-1; %#ok<PROPLC>
                Y = (1 : h) +y0-1; %#ok<PROPLC>
                [X, Y, binnedPic] = prepareSurfaceData(X, Y, pic);
                [obj.fitRes, obj.gof, obj.output] = fit([X, Y], binnedPic, f);
            end
            
            % set 2D params
            obj.OD = obj.fitRes.OD;
            obj.x0 = obj.fitRes.x0;
            obj.y0 = obj.fitRes.y0;
            obj.sigmaX = obj.fitRes.sigmaX;
            obj.sigmaY = obj.fitRes.sigmaY;
            obj.C = obj.fitRes.C;
            
%             % set light level into saveParamVal
%             calibratedAtom = (appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic + 1)*exp(-obj.C)-1; % background picture, calibrated to the atoms level using C parameter from the fit
%             appData.save.saveParamVal = mean(calibratedAtom(:));
%             set(appData.ui.etParamVal, 'String', appData.save.saveParamVal)
            
            %%%%%%%%%%% calculate residuals of the 2D fit:
%             if ( numel(pic) > 350^2 ) 
%                 residulas = obj.OD*( exp( -0.5 * (X - obj.x0).^2./obj.sigmaX.^2 - 0.5 * (Y - obj.y0).^2./obj.sigmaY.^2 ) ) + obj.C - binnedPic;
%             else
%                 residulas = obj.OD*( exp( -0.5 * (X - obj.x0).^2./obj.sigmaX.^2 - 0.5 * (Y - obj.y0).^2./obj.sigmaY.^2 ) ) + obj.C - pic;
%             end
%             disp(['Mean 2D Gaussian fit resiuals RMS = ' num2str(sqrt(mean(residulas(:).^2))) ])
%             figure;
%             imagesc(residulas)
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
           
            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
%            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
%                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
            [xData yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(obj.xCenter, obj.yCenter, binW);

            obj.xData = xData;
            obj.xStart = x0;
            obj.yData = yData;
            obj.yStart = y0;

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
           %            text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           text( 0, 0.8, 'fit function = A * {\ite}^{-(x-x_0)^2 / 2\sigma_x^2 - (y-y_0)^2 / 2\sigma_y^2} + C',...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.65,['A = ' num2str(obj.OD)], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.45, { ['x_0 = ' num2str(obj.x0 * appData.data.camera.xPixSz * 1000) ' mm;'], ...
               ['\sigma_x = ' num2str(obj.sigmaX * appData.data.camera.xPixSz * 1e6) ' um;']}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           
           text( 0.3, 0.45, {['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'],...
               ['\sigma_y = ' num2str(obj.sigmaY * appData.data.camera.yPixSz * 1e6) ' um']},...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.25,['C = ' num2str(obj.C)], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.15,['X/Y size aspect ratio = ' num2str(obj.sigmaX./obj.sigmaY)], 'fontsize', 12, 'Parent', appData.ui.pFitResults);
%            plot(rand(10,1));
%            title(['\fontsize{20} Look Big' char(10) ...
%                '\fontsize{10} Look small' char(10) ...
%                '\fontsize{20} Mixed_{\fontsize{8} underscore}'],'interpreter','tex');
           end
   end
end 
      
function ret = fitGauss2D_matrix( p, X, Y, g ) %#ok<DEFNU>
% amp = p(1);
% cx = p(2);
% cy = p(3);
% wx = p(4);
% wy = p(5);
% C = p(6)
ret = p(1)*( exp( -0.5 * (X - p(2)).^2./p(4).^2 - 0.5 * (Y - p(3)).^2./p(5).^2 ) ) + p(6) - g;
end
    
function ret = fitGauss2D_scalar( p, X, Y, g )
% amp = p(1);
% cx = p(2);
% cy = p(3);
% wx = p(4);
% wy = p(5);
% C = p(6)
ret = p(1)*( exp( -0.5 * (X - p(2)).^2./p(4).^2 - 0.5 * (Y - p(3)).^2./p(5).^2 ) ) + p(6) - g;
ret = sum(sum(ret.^2));
end