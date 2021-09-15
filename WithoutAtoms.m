classdef WithoutAtoms < PlotTypes
%WITHOUTATOMS Summary of this class goes here
%   Detailed explanation goes here

   properties
   end

   methods
       function appData = setPic(obj, appData, pic)
           obj.xStart = 1;
           obj.yStart = 1;
           obj.pic = pic;
           
           % last 
           appData.data.plots{appData.consts.plotTypes.withoutAtoms} = obj;
       end
             
       function [pic x0 y0] = getAnalysisPic(obj, appData)  %#ok<INUSD>
           [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.absorption}.getPic() ;
       end
       
       function [pic x0 y0] = getPic(obj) 
           x0 = obj.xStart;
           y0 = obj.yStart;
           pic = obj.pic;
       end
       
       function normalizedPic = normalizePic(obj, appData, pic)
           normalizedPic = pic/max(max(pic));
       end
       
       function [xData, yData] = getAnalysisXYDataVectors(obj, appData, x0, y0, avg) % return the vectors of the data
           [xData, yData] = appData.data.plots{appData.consts.plotTypes.absorption}.getXYDataVectors(x0, y0, avg);
       end
              
       function [xData, yData] = getXYDataVectors(obj, x0, y0, avg) % return the vectors of the data
           [h w] = size(obj.pic);
           xData = mean(obj.pic(max(y0-obj.yStart+1-avg, 1) : min(y0-obj.yStart+1+avg, h), :), 1);
           yData = mean(obj.pic(:, max(x0-obj.xStart+1-avg, 1) : min(x0-obj.xStart+1+avg, w)), 2)';
       end
       
   end
end 
