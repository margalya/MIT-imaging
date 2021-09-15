classdef PlotTypes
%PLOTTYPES Summary of this class goes here
%   Detailed explanation goes here

   properties ( SetAccess = protected )
       pic = [];
       xStart = -1;
       yStart = -1;
   end
   properties (Constant )
       lineColor = 0.85;
   end

   methods ( Abstract = true )
       appData = setPic(obj, appData, pic);
       
       [pic x0 y0] = getAnalysisPic(obj, appData) %return the correct pic for analysing according to the plot type
       [pic x0 y0] = getPic(obj) % return the correct pic according to the plot type
       
       normelizedPic = normalizePic(obj, appData, pic) % normelized to the maxVal
       [xData, yData] = getAnalysisXYDataVectors(obj, appData, x0, y0, avg) % return the vectors of the data
       [xData, yData] = getXYDataVectors(obj, x0, y0, avg) % return the vectors of the data
   end
   
   methods
       function squaredPic = createSquare(obj, appData, pic)
           squaredPic = pic;
       end
   end
end 
