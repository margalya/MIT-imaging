classdef FitXYCut < FitTypes
%FITONLYMAX Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'FitXYCut';
   end
   properties (SetAccess = private )
       stdv = 0;
   end

   methods 
       function appData = analyze(obj, appData) % do the analysis
           
           [x, y] = ginput(1);
           
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           [h w] = size(pic);
            xPosMax = round(x / (appData.data.camera.xPixSz * 1000)) - x0;
            yPosMax = round(y / (appData.data.camera.yPixSz * 1000)) - y0 + appData.data.camera.chipStart;
%            obj.maxVal = pic(yPosMax, xPosMax);
           obj.maxVal = max(max(pic));
           % center
           obj.xCenter = xPosMax + x0;
           obj.yCenter = yPosMax + y0;
           % unit size
           obj.xUnitSize = w;
           obj.yUnitSize = h;     
           
           obj.stdv = std(sqrt(pic(:).^2));
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           
           obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
           % last
           % set ROI pic - MUST be after defining ROI
           appData.data.fits{appData.consts.fitTypes.XYCut} = obj;
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
%            
%            % last 
%            appData.data.fits{appData.consts.fitTypes.onlyMaximum} = obj;
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x);
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           normalizedROI = pic(y, x);
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = pic/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = zeros(size(x));
           yFit = zeros(size(y));
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getPic();
           text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
           text( 50, 140, {['x_0 = ' num2str((obj.xCenter-x0+1) * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['y_0 = ' num2str((obj.yCenter-y0+1) * appData.data.camera.yPixSz * 1000) ' mm']}, ...
               'fontsize', 12);
           text( 50, 100, ['std = ' num2str(obj.stdv)], 'fontsize', 12); 
       end
       
   end
end 
