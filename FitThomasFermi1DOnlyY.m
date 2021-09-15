classdef FitThomasFermi1DOnlyY < FitTypes
%FITTHOMASFERMI1DONLYY Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'FitThomasFermi1DOnlyY';
   end
   properties (SetAccess = private )       
       y0 = -1;
       ODy = -1;
       TFhwY = -1; %Thomas Fermi half width
       Cy = -1;
       
       yRes = [];
       yGof = [];
       yOutput = [];
   end

   methods
       function appData = analyze(obj, appData) % do the analysis
           % 1D fit
           fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussianOnlyY};
           if ( isempty( fitObj.yRes ) )
               tmpFitType = appData.data.fitType;
               appData.data.fitType = appData.consts.fitTypes.oneDGaussianOnlyY;
               appData = appData.data.fits{appData.consts.fitTypes.oneDGaussianOnlyY}.analyze(appData);
               fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussianOnlyY};
               appData.data.fitType = tmpFitType;
           end
           
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);  % the whole pic
%            [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
           [h w] = size(pic);
           
           
           % 1D Thomas Fermi fit
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, fitObj.xCenter, fitObj.yCenter, appData.options.avgWidth);
            y = [1 : h] + y0-1;            
            
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', ...
                [1.5*fitObj.yRes.Ay fitObj.yRes.y0 fitObj.yRes.sigmaY fitObj.yRes.Cy] , ...
                'lower', [0 0 0  -inf]);
            f = fittype('ODy * max( 1 - ( (y-y0)./TFhwY ).^2, 0 ).^(3/2) + Cy', ...
                'coefficients', {'ODy', 'y0', 'TFhwY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
            
            
            % set 1D params
           obj.y0 = obj.yRes.y0;

           obj.ODy = obj.yRes.ODy;
           obj.TFhwY = obj.yRes.TFhwY; %Thomas Fermi half width
           obj.Cy = obj.yRes.Cy;
            
            %set fit params
            obj.xCenter = round(fitObj.xCenter); % should be indexes (integers)
            obj.yCenter= round(obj.y0);
            obj.xUnitSize = w;
            obj.yUnitSize = obj.TFhwY;
            obj.maxVal = obj.ODy;
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData);
%            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           
           obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
%             [xData yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(obj.xCenter, obj.yCenter, appData.options.avgWidth);
            
            obj.xData = xData;
            obj.xStart = x0;
            obj.yData = yData;
            obj.yStart = y0;
            
           % last 
           appData.data.fits{appData.consts.fitTypes.oneDTFOnlyY} = obj;
%            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - obj.Cy;
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
%            [X, Y] = meshgrid(x, y);
% %            normalizedROI = obj.OD * max(1 - ( (X-obj.x0)./ (2*obj.TFhwX) ).^2 - ( (Y-obj.y0)./(2*obj.TFhwY) ).^2 ,0).^(3/2);
%            normalizedROI =obj.ODx * max(1  - ( (Y-obj.y0)./(obj.TFhwY) ).^2 ,0).^(3/2);
            normalizedROI = pic(y, x);
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic - obj.Cy) / obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = zeros(size(x));
           yFit = obj.ODy * max( 1 - ( (y-obj.y0)./(obj.TFhwY) ).^2, 0 ).^(3/2) + obj.Cy;
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           chipStart = appData.data.camera.chipStart;
           cla(appData.ui.pFitResults)
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           text( 0, 0.8, {'fit func = A_y * max(1 - (y-y_0)^2/w_y^2, 0)^{3/2}+ Cy'}, 'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.4, {['A_y = ' num2str(obj.ODy)], ...
               ['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...                
               ['w_y = ' num2str(obj.TFhwY * appData.data.camera.yPixSz * 1000) ' mm'], ...
               ['C_y = ' num2str(obj.Cy)]}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults); 
       end
   end
end 
  
