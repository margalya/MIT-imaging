classdef Size < ROIUnits
% MM Summary of this class goes here
%   Detailed explanation goes here

   properties
   end
% *appData.data.camera.xPixSz * 1000
   methods 
       function [x0 y0 x1 y1] = getROICoords(obj, appData, fitObj)
           xSize = appData.data.ROISizeX / 1000 / appData.data.camera.xPixSz;
           ySize = appData.data.ROISizeY / 1000 / appData.data.camera.yPixSz;
%            [pic xStart yStart] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData); 
           [pic xStart yStart] = appData.data.plots{appData.consts.plotTypes.absorption}.getPic();
           [h w] = size(pic);
           
           xCenter = appData.data.ROICenterX / 1000 / appData.data.camera.xPixSz + xStart-1; % in pixels
           yCenter = appData.data.ROICenterY / 1000 / appData.data.camera.yPixSz + yStart-1;           
           xCenter = max(round(xCenter), xStart);
           yCenter = max(round(yCenter), yStart);
           xCenter = min(round(xCenter), xStart+w-1);
           yCenter = min(round(yCenter), yStart+h-1);
           
           x0 = xCenter - xSize;
           x1 = xCenter + xSize;
           y0 = yCenter - ySize;
           y1 = yCenter + ySize;           
           
           % this if...else... is used when no fit was done - dosn't exist
           % in ROI Size (when the ROI is already known)
%            if ( fitObj.xCenter == -1 || fitObj.yCenter == -1 )
%                x0 = xStart;
%                x1 = xStart+w-1;
%                y0 = yStart;
%                y1 = yStart+h-1;  
%            else
               x0 = max(round(x0), xStart);
               x1 = min(round(x1), xStart+w-1);
               y0 = max(round(y0), yStart);
               y1 = min(round(y1), yStart+h-1);
%            end
       end
       
       function [x0 y0] = getCenter(obj, appData, fitObj) %#ok<INUSD>
           xCenter = appData.data.ROICenterX / 1000 / appData.data.camera.xPixSz; % in pixels
           yCenter = appData.data.ROICenterY / 1000 / appData.data.camera.yPixSz;           
           [pic xStart yStart] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData); 
           [h w] = size(pic);
           xCenter = max(round(xCenter), xStart);
           yCenter = max(round(yCenter), yStart);
           xCenter = min(round(xCenter), xStart+w-1);
           yCenter = min(round(yCenter), yStart+h-1);
           x0 = xCenter; %TODO: round, min max
           y0 = yCenter;
           x0 = fitObj.xCenter; %TODO: round, min max
           y0 = fitObj.yCenter;
       end
   end
end 
