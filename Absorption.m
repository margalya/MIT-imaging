classdef Absorption < PlotTypes
%ABSORPTION Summary of this class goes here
%   Detailed explanation goes here

   properties
   end

   methods
       function appData = setPic(obj, appData, pic)
           obj.xStart = 1;
           obj.yStart = appData.data.camera.chipStart;
           obj.pic = pic(obj.yStart : end, :);
           
           % last 
           appData.data.plots{appData.consts.plotTypes.absorption} = obj;
       end
             
       function [pic x0 y0] = getAnalysisPic(obj, appData)  %#ok<INUSD>
           [pic x0 y0] = obj.getPic() ;
       end
       
       function [pic x0 y0] = getPic(obj) 
           x0 = obj.xStart;
           y0 = obj.yStart;
           pic = obj.pic;
%            [y1 x1] = size(obj.pic);
%            pic = obj.pic([0:y1]+y0, [0:x1]+x0);
       end
       
       function normalizedPic = normalizePic(obj, appData, pic)
           normalizedPic = appData.data.fits{appData.data.fitType}.normalizePic(pic);
       end
       
       function [xData, yData] = getAnalysisXYDataVectors(obj, appData, x0, y0, avg) % return the vectors of the data
           [xData, yData] = obj.getXYDataVectors(x0, y0, avg);
       end
       
       function [xData, yData] = getXYDataVectors(obj, x0, y0, avg) % return the vectors of the data
           [h w] = size(obj.pic);
%            xData = mean(obj.pic(min(y0-avg, 1) : max(y0+avg, h), :), 1);
%            yData = mean(obj.pic(:, min(x0-avg, 1) : max(x0+avg, w)), 2)';

           xData = mean(obj.pic(max(y0-obj.yStart+1-avg, 1) : min(y0-obj.yStart+1+avg, h), :), 1);
           yData = mean(obj.pic(:, max(x0-obj.xStart+1-avg, 1) : min(x0-obj.xStart+1+avg, w)), 2)';
           
%                 yData = mean(appData.data.pic(appData.data.yPosMax-appData.data.ROIBottom-avg : ...
%                     appData.data.yPosMax-appData.data.ROIBottom+avg, x), 1);
%                 xData = mean(appData.data.pic( y, appData.data.xPosMax-appData.data.ROILeft-avg : ...
%                     appData.data.xPosMax-appData.data.ROILeft+avg), 2)'; 
       end
       
       function squaredPic = createSquare(obj, appData, pic)
            fitObj = appData.data.fits{appData.data.fitType};
            squaredPic = pic;
            squaredPic([fitObj.ROIBottom-1 : fitObj.ROIBottom] -obj.yStart+1, [fitObj.ROILeft : fitObj.ROIRight] -obj.xStart+1) = ...
                ones(2, fitObj.ROIRight - fitObj.ROILeft+1) * obj.lineColor;
            squaredPic([fitObj.ROITop : fitObj.ROITop+1] -obj.yStart+1, [fitObj.ROILeft : fitObj.ROIRight] -obj.xStart+1) = ...
                ones(2, fitObj.ROIRight - fitObj.ROILeft+1) * obj.lineColor;
            squaredPic([fitObj.ROITop : fitObj.ROIBottom] -obj.yStart+1, [fitObj.ROILeft : fitObj.ROILeft+1] -obj.xStart+1) = ...
                ones(fitObj.ROIBottom - fitObj.ROITop+1, 2) * obj.lineColor;
            squaredPic([fitObj.ROITop : fitObj.ROIBottom] -obj.yStart+1, [fitObj.ROIRight-1 : fitObj.ROIRight] -obj.xStart+1) = ...
                ones(fitObj.ROIBottom - fitObj.ROITop+1, 2) * obj.lineColor;
       end
       
%        function [xData, yData] = getXYData(obj, appData)
%            pic = appData.data.absorption;
%            xMax = appData.data.fits{appData.data.fitType}.xPosMax;
%            yMax = appData.data.fits{appData.data.fitType}.yPosMax;  
%            avg = appData.options.avgWidth;
%            
%            xData = mean(pic(yMax-avg : yMax+avg, :), 2)';
%            yData = mean(pic(:, xMax-avg : xMax+avg), 1);           
%        end
   end
end 
