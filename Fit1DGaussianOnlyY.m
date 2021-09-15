classdef Fit1DGaussianOnlyY < FitTypes
%FIT1DGAUSSIANONLYY Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'Fit1DGaussianOnlyY';
   end
   properties (SetAccess = private )
%        xRes = [];
%        xGof = [];
%        xOutput = [];
       yRes = [];
       yGof = [];
       yOutput = [];
   end

   methods 
       function appData = analyze(obj, appData) % do the analysis
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           [h w] = size(pic);
           %binning and max
           binnedData = binning ( pic, appData.options.avgWidth*2+1);
           [maxes, indexes] = max(binnedData);                     % find maximum
           [maxValue, xPosMax] = max(maxes);
           yPosMax = indexes(xPosMax);   
%            obj.maxVal = maxValue;% / (appData.options.avgWidth*2)^2; %no need, already done in binning.m
           % center
           xCenter = (appData.options.avgWidth*2+1) * (xPosMax ) + x0-appData.options.avgWidth;
           yCenter = (appData.options.avgWidth*2+1) * (yPosMax ) + y0-appData.options.avgWidth;
%            obj.xCenter = appData.options.avgWidth*2 * (xPosMax - 0.5) + x0;
%            obj.yCenter = appData.options.avgWidth*2 * (yPosMax - 0.5) + y0; 
%            % unit size
%            obj.xUnitSize = w/2;
%            obj.yUnitSize = h/2;   

            % fitting
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
%             x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;            
%             [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', 'gauss1'); 
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', 'gauss1');  
%             xCenter = round(obj.xRes.b1);
            yCenter = round(obj.yRes.b1);
            
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
%             s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes.a1 obj.xRes.b1 obj.xRes.c1/sqrt(2)]);%[1 width/2 width/10]);
%             f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))', 'coefficients', {'Ax', 'x0', 'sigmaX'}, 'independent', 'x', 'dependent', 'y', 'options', s);
%             [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes.a1 obj.yRes.b1 obj.yRes.c1/sqrt(2)]);%[1 width/2 width/10]);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))', 'coefficients', {'Ay', 'y0', 'sigmaY'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
%             xCenter = round(obj.xRes.x0);
            yCenter = round(obj.yRes.y0);
%              yCenter = obj.yRes.y0;
             
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
%             s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes.Ax obj.xRes.x0 obj.xRes.sigmaX 0], 'Lower', [0 -Inf 0 -Inf]);
%             f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
%             [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes.Ay obj.yRes.y0 obj.yRes.sigmaY 0], 'Lower', [0 -Inf 0 -Inf]);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 

            
%             i = 0;
%             while ( ( (abs(round(obj.xRes.x0)-xCenter) > 1) || (abs(round(obj.yRes.y0)-yCenter) > 1) ) && i < 10 )
%                 i = i + 1;
%                 xCenter = round(obj.xRes.x0);
%                 yCenter = round(obj.yRes.y0);
% %                 appData.data.xPosMax = round(xRes.x0);
% %                 appData.data.yPosMax = round(yRes.y0);
% 
%                 [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
%                 s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes.Ax obj.xRes.x0 obj.xRes.sigmaX obj.xRes.Cx]);
%                 f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
%                 [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
%                 s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes.Ay obj.yRes.y0 obj.yRes.sigmaY obj.yRes.Cy]);
%                 f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
%                 [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
%             end
            
            %set fit params
            obj.xCenter = round(xCenter); % should be indexes (integers)
            obj.yCenter= obj.yRes.y0;
            obj.xUnitSize = w;
            obj.yUnitSize = obj.yRes.sigmaY;
            obj.maxVal = obj.yRes.Ay;
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           
           obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
           % last 
           appData.data.fits{appData.consts.fitTypes.oneDGaussianOnlyY} = obj;
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - obj.yRes.Cy;
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y) %#ok<INUSL>
%            [X, Y] = meshgrid(x, y);
%            normalizedROI = mean([obj.xRes.Ax obj.yRes.Ay]) * exp( -0.5 * ((X-obj.xRes.x0).^2 / obj.xRes.sigmaX^2 + (Y-obj.yRes.y0).^2 / obj.yRes.sigmaY^2) );
            normalizedROI = pic(y, x);
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic - obj.yRes.Cy)/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = zeros(size(x));
           yFit = obj.yRes.Ay * exp( -0.5 * (y-obj.yRes.y0).^2 / obj.yRes.sigmaY^2 ) + obj.yRes.Cy;
       end
       
       function  plotFitResults(obj, appData)  % plots the text           
           chipStart = appData.data.camera.chipStart;
           cla(appData.ui.pFitResults)
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           text( 0, 0.8, 'fit function = A_y * {\ite}^{-(y-y_0)^2 / 2\sigma_y^2}+C_y', 'fontsize', 12, 'Parent', appData.ui.pFitResults);
%            text( 50, 80, {['A_x = ' num2str(obj.xRes.Ax)], ...
%                ['x_0 = ' num2str((obj.xRes.x0) * appData.data.camera.xPixSz * 1000) ' mm'], ...
%                ['\sigma_x = ' num2str(obj.xRes.sigmaX * appData.data.camera.xPixSz * 1000) ' mm'], ...
%                ['C_x = ' num2str(obj.xRes.Cx)], ...
%                ['R^2 = ' num2str(obj.xGof.rsquare)]}, ...
%                'fontsize', 12);
           text( 0, 0.4, {['A_y = ' num2str(obj.yRes.Ay)], ...
               ['y_0 = ' num2str((obj.yRes.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...
               ['\sigma_y = ' num2str(obj.yRes.sigmaY * appData.data.camera.yPixSz * 1000) ' mm'], ...
               ['C_y = ' num2str(obj.yRes.Cy)], ...
               ['R^2 = ' num2str(obj.yGof.rsquare)]}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
       end
       
   end
end 
