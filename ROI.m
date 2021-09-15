classdef ROI < PlotTypes
%ROI Summary of this class goes here
%   Detailed explanation goes here

   properties
   end

   methods
       function appData = setPic(obj, appData, pic)
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.absorption}.getPic();
           fitObj = appData.data.fits{appData.data.fitType};
           
           [ROILeft ROITop ROIRight ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, fitObj);
           obj.xStart = ROILeft;
           obj.yStart = ROITop;
           obj.pic = pic( [obj.yStart : ROIBottom] -y0+1, [obj.xStart : ROIRight] -x0+1);
           
           if ( obj.xStart == -1 && obj.yStart == -1 && isempty(obj.pic) )
               obj.xStart = x0;
               obj.yStart = y0;
               obj.pic = pic;
           end               
           
           % last 
           appData.data.plots{appData.consts.plotTypes.ROI} = obj;
       end
             
       function [pic x0 y0] = getAnalysisPic(obj, appData)  %#ok<INUSD>
%            if ( isempty(obj.pic) == 1 )
%                [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.absorption}.getPic();
%            else
               [pic x0 y0] = obj.getPic();
%            end
       end
       
       function [pic x0 y0] = getPic(obj) 
%            if ( isempty(obj.pic) == 1 )
%                [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.absorption}.getPic();
%            else
               x0 = obj.xStart;
               y0 = obj.yStart;
               pic = obj.pic;
%            end
       end
       
       function normalizedPic = normalizePic(obj, appData, pic)
           normalizedPic = appData.data.fits{appData.data.fitType}.normalizePic(pic);
       end
       
       function [xData, yData] = getAnalysisXYDataVectors(obj, appData, x0, y0, avg) % return the vectors of the data
           [xData, yData] = obj.getXYDataVectors(x0, y0, avg);
       end
       
       function [xData, yData] = getXYDataVectors(obj, x0, y0, avg) % return the vectors of the data
           [h w] = size(obj.pic);
           xData = mean(obj.pic(max(y0-obj.yStart+1-avg, 1) : min(y0-obj.yStart+1+avg, h), :), 1);
           yData = mean(obj.pic(:, max(x0-obj.xStart+1-avg, 1) : min(x0-obj.xStart+1+avg, w)), 2)';
       end
      
   end
end 
