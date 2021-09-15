classdef FitOnlyMax < FitTypes
%FITONLYMAX Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'FitOnlyMax';
   end
   properties ( SetAccess = private)
       stdv = 0;
       meanCountWith = -1; % mean photons count of with atoms image (in ROI)
       meanCountWithout = -1; % mean photon count without (in ROI)
       meanCountDark = -1;
       countRatio = -1; % ratio of counts
       photonCount = []; % number of collected photons, fluorescence imaging
   end

   methods 
       function appData = analyze(obj, appData) % do the analysis
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           [h w] = size(pic);
           %binning and max
%            binW = appData.options.avgWidth;
           %            binnedData = binning ( pic, binW*2+1);
           binnedData = LowPassFilter(pic,30,appData.options.avgWidth);
           [maxes, indexes] = max(binnedData);                     % find maximum
           [maxValue, xPosMax] = max(maxes);
           yPosMax = indexes(xPosMax);
%            obj.maxVal = max(pic(:));
           obj.maxVal = maxValue;% / (appData.options.avgWidth*2)^2; %no need, already done in binning.m
           % center
           obj.xCenter = xPosMax + x0; %(binW*2+1) * (xPosMax ) + x0-binW;
           obj.yCenter = yPosMax + y0; %(binW*2+1) * (yPosMax ) + y0-binW; 
           % unit size
           obj.xUnitSize = w;
           obj.yUnitSize = h;   
           
%            obj.stdv = std(sqrt(pic(:).^2));
           obj.stdv = std(pic(:));
           obj.meanCountWith = mean2(appData.data.plots{appData.consts.plotTypes.withAtoms}.pic);
           obj.meanCountWithout = mean2(appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic);
           obj.countRatio = obj.meanCountWith ./ obj.meanCountWithout;
           obj.meanCountDark = mean2(appData.data.plots{appData.consts.plotTypes.dark}.pic);
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           
           obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
           if isfield(appData.consts.cameras{appData.options.cameraType}, 'photonPerADU')
               obj.photonCount = round(sum(pic(:)) * appData.consts.cameras{appData.options.cameraType}.photonPerADU); %photon counts
           end
           
           [xData yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(obj.xCenter, obj.yCenter, 0);
           
           obj.xData = xData;
           obj.xStart = x0;
           obj.yData = yData;
           obj.yStart = y0;
           
           % last
           % set ROI pic - MUST be after defining ROI
           appData.data.fits{appData.consts.fitTypes.onlyMaximum} = obj;
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
%            picMean = mean(pic(:));
%            picStd = std(pic(:));
%            pic = (pic - picMean - 2*picStd);
           normalizedPic = pic/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = zeros(size(x));
           yFit = zeros(size(y));
       end
       
       function  plotFitResults(obj, appData)  % plots the text
           [pic, ~, ~] = appData.data.plots{appData.data.plotType}.getPic();
%            set(appData.ui.pFitResults, 'String', ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 18);
           cla(appData.ui.pFitResults)
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 18, 'Parent', appData.ui.pFitResults);
           if ~isempty(obj.photonCount)
               text( 0, 0.85, ['Photon count: ' num2str(obj.photonCount)], 'fontSize', 12, 'Parent', appData.ui.pFitResults);
               photonPerPixel = obj.photonCount./numel(pic);
               text( 0, 0.75, ['Mean photons per pixel = ' num2str(photonPerPixel) ], 'fontSize', 12, 'Parent', appData.ui.pFitResults);
           end
           meanWith = round( mean(appData.data.plots{appData.consts.plotTypes.withAtoms}.pic(:)) * 10) ./ 10 ;
           meanWithout = round( mean(appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic(:)) * 10) ./ 10;
           meanDark = round( mean(appData.data.plots{appData.consts.plotTypes.dark}.pic(:)) * 10) ./ 10;
           maxWithout = max(appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic(:));
           
           textToDisp = {['OD_ = ' num2str(obj.maxVal) ], ...
               ['x_0 = ' num2str((obj.xCenter) * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['y_0 = ' num2str((obj.yCenter-appData.data.camera.chipStart) * appData.data.camera.yPixSz * 1000) ' mm'],...
               ['std = ' num2str(obj.stdv)],...
               ['Mean image OD = ' num2str(mean(pic(:)))],...
               ['Mean count with, w/o, dark = ' num2str(meanWith) ', ' num2str(meanWithout) ', ' num2str(meanDark)],...
               ['Max pixel count (w/o) = ' num2str(maxWithout) ]};
           
           text( 0, 0.35, textToDisp, 'fontsize', 11, 'Parent', appData.ui.pFitResults);
           
       end
       
   end
end 
