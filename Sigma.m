classdef Sigma < ROIUnits
% SIGMA Summary of this class goes here
%   Detailed explanation goes here

   properties
   end

   methods 
       function [x0 y0 x1 y1] = getROICoords(obj, appData, fitObj)
           xSize = appData.data.ROISizeX;
           ySize = appData.data.ROISizeY;
           x0 = fitObj.xCenter - xSize*fitObj.xUnitSize;
           x1 = fitObj.xCenter + xSize*fitObj.xUnitSize;
           y0 = fitObj.yCenter - ySize*fitObj.yUnitSize;
           y1 = fitObj.yCenter + ySize*fitObj.yUnitSize;           
           
%            [pic xStart yStart] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData); 
           [pic xStart yStart] = appData.data.plots{appData.consts.plotTypes.absorption}.getAnalysisPic(appData);
           [h w] = size(pic);
           
           % this if...else... is used when no fit was done - dosn't exist
           % in ROI Size (when the ROI is already known)
           if ( fitObj.xUnitSize == -1 || fitObj.yUnitSize == -1 )
               x0 = xStart;
               x1 = xStart+w-1;
               y0 = yStart;
               y1 = yStart+h-1;  
           else
               x0 = max(round(x0), xStart);
               x1 = min(round(x1), xStart+w-1);
               y0 = max(round(y0), yStart);
               y1 = min(round(y1), yStart+h-1);
           end
       end
       
       function [x0 y0] = getCenter(obj, appData, fitObj)
           x0 = fitObj.xCenter;%TODO: round, min max
           y0 = fitObj.yCenter;
       end
   end
end 
