classdef FitBiModal1D < FitTypes
%FITBIMODAL1D Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'FitBiModal1D';
   end
   properties (SetAccess = private )
       x0 = -1;
       y0 = -1;
       
       ampTFx = -1;
       TFhwX = -1; %Thomas Fermi half width
       ampGx = -1;
       sigmaX = -1; %Gaussian sigma
       Cx = -1;
       ampTFy = -1;
       TFhwY = -1; %Thomas Fermi half width
       ampGy = -1;
       sigmaY = -1; %Gaussian sigma
       Cy = -1;
       
       xRes = [];
       xGof = [];
       xOutput = [];
       yRes = [];
       yGof = [];
       yOutput = [];
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
           
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);  % the whole pic
%            [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
           [h w] = size(pic);
           
           
           % 1D BiModal fit
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, fitObj.xCenter, fitObj.yCenter, appData.options.avgWidth);
            x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;            
            
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', ...
                [1*fitObj.xRes.Ax fitObj.xRes.x0 0.5*fitObj.xRes.sigmaX 1*fitObj.xRes.Ax 10*fitObj.xRes.sigmaX, fitObj.xRes.Cx], ...
                'lower', [0 -inf 0 0 0 -inf]);
            f = fittype('ampTFx * max( 1 - ( (x-x0)./TFhwX ).^2, 0 ).^(3/2) + ampGx*exp( -0.5 * ( (x-x0).^2 / sigmaX^2 ) ) + Cx', ...
                'coefficients', {'ampTFx', 'x0', 'TFhwX', 'ampGx', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
            
%             s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', ...
%                 [1.5*fitObj.yRes.Ay fitObj.yRes.y0 0.25*fitObj.yRes.sigmaY 0.1*fitObj.yRes.Ay 10*fitObj.yRes.sigmaY, fitObj.yRes.Cy] , ...
%                 'lower', [0 -inf 0 0 0 -inf]);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', ...
                [obj.xRes.ampTFx fitObj.yRes.y0 obj.xRes.TFhwX obj.xRes.ampGx obj.xRes.sigmaX, fitObj.yRes.Cy] , ...
                'lower', [0 -inf 0 0 0 -inf]);
            
            f = fittype('ampTFy * max( 1 - ( (y-y0)./TFhwY ).^2, 0 ).^(3/2) + ampGy*exp( -0.5 * ( (y-y0).^2 / sigmaY^2 ) ) + Cy', ...
                'coefficients', {'ampTFy', 'y0', 'TFhwY', 'ampGy', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
            
            
            % set 1D params
           obj.x0 = obj.xRes.x0;
           obj.y0 = obj.yRes.y0;

           obj.ampTFx = obj.xRes.ampTFx;
           obj.TFhwX = obj.xRes.TFhwX; %Thomas Fermi half width
           obj.ampGx = obj.xRes.ampGx;
           obj.sigmaX = obj.xRes.sigmaX; %Gaussian sigma
           obj.Cx = obj.xRes.Cx;
           obj.ampTFy = obj.yRes.ampTFy;
           obj.TFhwY = obj.yRes.TFhwY; %Thomas Fermi half width
           obj.ampGy = obj.yRes.ampGy;
           obj.sigmaY = obj.yRes.sigmaY; %Gaussian sigma
           obj.Cy = obj.yRes.Cy;
            
            %set fit params
            obj.xCenter = round(obj.x0); % should be indexes (integers)
            obj.yCenter= round(obj.y0);
%             obj.xUnitSize = obj.TFhwX;
%             obj.yUnitSize = obj.TFhwY;
            obj.xUnitSize = obj.sigmaX;
            obj.yUnitSize = obj.sigmaY;
            obj.maxVal = min([obj.ampTFx+obj.ampGx obj.ampTFy+obj.ampGy]);%mean([obj.ampTFx obj.ampTFy]) + mean([obj.ampGx obj.ampGy]);
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData);
%            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           
           obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
           obj.xData = xData;
           obj.xStart = x0;
           obj.yData = yData;
           obj.yStart = y0;
           
           % last 
           appData.data.fits{appData.consts.fitTypes.oneDBiModal} = obj;
%            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - mean([obj.xRes.Cx obj.yRes.Cy]);
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
%            normalizedROI = obj.OD * max(1 - ( (X-obj.x0)./ (2*obj.TFhwX) ).^2 - ( (Y-obj.y0)./(2*obj.TFhwY) ).^2 ,0).^(3/2);
           normalizedROI = mean([obj.ampTFx obj.ampTFy]) * max(1 - ( (X-obj.x0)./ (obj.TFhwX) ).^2 - ( (Y-obj.y0)./(obj.TFhwY) ).^2 ,0).^(3/2) + ...
                mean([obj.ampGx obj.ampGy])*exp( -0.5 * ((X-obj.x0).^2 / obj.sigmaX^2 + (Y-obj.y0).^2 / obj.sigmaY^2) );
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-mean([obj.Cx obj.Cy]))/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit1 = obj.ampGx*exp( -0.5 * ( (x-obj.x0).^2 / obj.sigmaX^2 ) ) + obj.Cx;
           xFit2 = obj.ampTFx * max( 1 - ( (x-obj.x0)./(obj.TFhwX) ).^2, 0 ).^(3/2) + xFit1;
           yFit1 = obj.ampGy*exp( -0.5 * ( (y-obj.y0).^2 / obj.sigmaY^2) ) + obj.Cy;
           yFit2 = obj.ampTFy * max( 1 - ( (y-obj.y0)./(obj.TFhwY) ).^2, 0 ).^(3/2) + yFit1;
           
           xFit = [xFit2; xFit1];
           yFit = [yFit2; yFit1];
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           chipStart = appData.data.camera.chipStart;
           cla(appData.ui.pFitResults)
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 16, 'Parent', appData.ui.pFitResults);%/1e6) '*10^6'], 'fontSize', 20);
           text( 0, 0.8, {'fit func = A_{TFx} * max(1 - (x-x_0)^2/w_x^2, 0)^{3/2}+  A_{Gx}*{\ite}^{ -(x-x_0)^2 / 2\sigma_x^2} + Cx' ... 
               '            + A_{TFy} * max(1 - (y-y_0)^2/w_y^2, 0)^{3/2}+  A_{Gy}*{\ite}^{ -(y-y_0)^2 / 2\sigma_y^2} + Cy'}, 'fontsize', 10, 'Parent', appData.ui.pFitResults);
%            text( 50, 110, ['A_{TFx}/A_{Gx} = ' num2str(obj.ampTFx) ' / ' num2str(obj.ampGx) ], 'fontsize', 12); 
           text( 0, 0.35, {['A_{TFx}/A_{Gx} = ' num2str(obj.ampTFx,'%.3f') ' / ' num2str(obj.ampGx,'%.3f') ], ...
               ['x_0 = ' num2str(obj.x0 * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['w_x = ' num2str(obj.TFhwX * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['\sigma_x = ' num2str(obj.sigmaX * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['C_x = ' num2str(obj.Cx)]}, ...
               'fontsize', 10, 'Parent', appData.ui.pFitResults);
           text( 0.4, 0.35, {['A_{TFy}/A_{Gy} = ' num2str(obj.ampTFy,'%.3f') ' / ' num2str(obj.ampGy,'%.3f') ], ...
               ['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...                
               ['w_y = ' num2str(obj.TFhwY * appData.data.camera.yPixSz * 1000) ' mm'], ...
               ['\sigma_y = ' num2str(obj.sigmaY * appData.data.camera.yPixSz * 1000) ' mm'], ...
               ['C_y = ' num2str(obj.Cy)]}, ...
               'fontsize', 10, 'Parent', appData.ui.pFitResults);
%            text( 0, 0.3, ['Mean amp = ' num2str((obj.ampTFx + obj.ampTFy + obj.ampGx + obj.ampGy)./2)], 'fontsize', 10, 'Parent', appData.ui.pFitResults); 
       end
   end
end 
  
