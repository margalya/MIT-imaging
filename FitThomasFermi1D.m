classdef FitThomasFermi1D < FitTypes
%FITTHOMASFERMI1D Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'FitThomasFermi1D';
   end
   properties (SetAccess = private )
       x0 = -1;
       ODx = -1;
       TFhwX = -1; %Thomas Fermi half width
       Cx = -1;
       
       y0 = -1;
       ODy = -1;
       TFhwY = -1; %Thomas Fermi half width
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
           
           
           % 1D Thomas Fermi fit
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, fitObj.xCenter, fitObj.yCenter, appData.options.avgWidth);
            x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;            
            
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', ...
                [fitObj.xRes.Ax fitObj.xRes.x0 fitObj.xRes.sigmaX  fitObj.xRes.Cx], ...
                'lower', [0 0 0 -inf]);
            f = fittype('ODx * max( 1 - ( (x-x0)./TFhwX ).^2, 0 ).^(3/2)  + Cx', ...
                'coefficients', {'ODx', 'x0', 'TFhwX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
            
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', ...
                [1.5*fitObj.yRes.Ay fitObj.yRes.y0 fitObj.yRes.sigmaY fitObj.yRes.Cy] , ...
                'lower', [0 0 0  -inf]);
            f = fittype('ODy * max( 1 - ( (y-y0)./TFhwY ).^2, 0 ).^(3/2) + Cy', ...
                'coefficients', {'ODy', 'y0', 'TFhwY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
            
            
            % set 1D params
           obj.x0 = obj.xRes.x0;
           obj.y0 = obj.yRes.y0;

           obj.ODx = obj.xRes.ODx;
           obj.TFhwX = obj.xRes.TFhwX; %Thomas Fermi half width
           obj.Cx = obj.xRes.Cx;
           obj.ODy = obj.yRes.ODy;
           obj.TFhwY = obj.yRes.TFhwY; %Thomas Fermi half width
           obj.Cy = obj.yRes.Cy;
            
            %set fit params
            obj.xCenter = round(obj.x0); % should be indexes (integers)
            obj.yCenter= round(obj.y0);
            obj.xUnitSize = obj.TFhwX;
            obj.yUnitSize = obj.TFhwY;
            obj.maxVal = mean([obj.ODx obj.ODy]);
           
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
           appData.data.fits{appData.consts.fitTypes.oneDTF} = obj;
%            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - mean([obj.xRes.Cx obj.yRes.Cy]);
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
%            normalizedROI = obj.OD * max(1 - ( (X-obj.x0)./ (2*obj.TFhwX) ).^2 - ( (Y-obj.y0)./(2*obj.TFhwY) ).^2 ,0).^(3/2);
           normalizedROI = mean([obj.ODx obj.ODy]) * max(1 - ( (X-obj.x0)./ (obj.TFhwX) ).^2 - ( (Y-obj.y0)./(obj.TFhwY) ).^2 ,0).^(3/2) ;
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-mean([obj.Cx obj.Cy]))/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = obj.ODx * max( 1 - ( (x-obj.x0)./(obj.TFhwX) ).^2, 0 ).^(3/2) + obj.Cx;
           yFit = obj.ODy * max( 1 - ( (y-obj.y0)./(obj.TFhwY) ).^2, 0 ).^(3/2) + obj.Cy;
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           chipStart = appData.data.camera.chipStart;
           cla(appData.ui.pFitResults)
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           text( 0, 0.75, {'fit func = A_x * max(1 - (x-x_0)^2/w_x^2, 0)^{3/2}+ Cx' ... 
               '            + A_y * max(1 - (y-y_0)^2/w_y^2, 0)^{3/2}+ Cy'}, 'fontsize', 12, 'Parent', appData.ui.pFitResults);
%            text( 50, 110, ['A_{TFx}/A_{Gx} = ' num2str(obj.ampTFx) ' / ' num2str(obj.ampGx) ], 'fontsize', 12); 
           text( 0, 0.3, {['A_x = ' num2str(obj.ODx) ], ...
               ['x_0 = ' num2str(obj.x0 * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['w_x = ' num2str(obj.TFhwX * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['C_x = ' num2str(obj.Cx)]}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0.3, 0.3, {['A_y = ' num2str(obj.ODy)], ...
               ['y_0 = ' num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...                
               ['w_y = ' num2str(obj.TFhwY * appData.data.camera.yPixSz * 1000) ' mm'], ...
               ['C_y = ' num2str(obj.Cy)]}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
%            text( 50, 20, ['C = ' num2str(obj.C)], 'fontsize', 12); 
       end
   end
end 
  
